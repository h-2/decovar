#include <cstddef>
#include <sharg/all.hpp>
#include <variant>

#include <bio/io/exception.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/reader.hpp>
#include <bio/io/var/writer.hpp>

using record_t = bio::io::var::record_default;
using header_t = bio::io::var::header;

struct decovar_error : std::runtime_error
{
    decovar_error(std::string msg) : std::runtime_error{msg} {}
};

/* vectors for fields of multiplicity A, R or G indicating whether that value at the positions shall be removed(1)
* or not (0) */
struct filter_vectors_t
{
    std::vector<int /*bool*/> A;
    std::vector<int /*bool*/> R;
    std::vector<int /*bool*/> G;
};

struct program_options
{
    std::string input_file;
    std::string output_file = "-";

    float rare_af_threshold = 1e-5;
};

program_options parse_options(int const argc, char const * const * const argv)
{
    program_options o;

    sharg::parser parser{"deCoVar", argc, argv};
    parser.info.author            = "Hannes Hauswedell";
    parser.info.short_description = "Reduce allele complexity in VCF files.";
    parser.info.version           = "0.1.0";
    parser.add_positional_option(o.input_file,
                                 sharg::config{.description = "Path to input file or '-' for stdin.",
                                               .required    = true,
                                               .validator   = sharg::input_file_validator{{"vcf", "vcf.gz", "bcf"}}});
    parser.add_option(o.output_file,
                      sharg::config{
                        .short_id    = 'o',
                        .long_id     = "output",
                        .description = "Path to output file or '-' for stdout.",
                        .validator   = sharg::output_file_validator{sharg::output_file_open_options::create_new,
                                                                    {"vcf", "vcf.gz", "bcf"}}
    });

    parser.add_option(o.rare_af_threshold,
                      sharg::config{
                        .long_id     = "rare-af-thresh",
                        .description = "Threshold below which alleles are considered rare.",
                        .validator   = sharg::arithmetic_range_validator{1e-10, 1e-1}
    });

    parser.parse();
    return o;
}

void determine_filtered_alleles(record_t::info_t const & record_info,
                                size_t const             record_no,
                                size_t const             n_alts,
                                program_options const &  opts,
                                filter_vectors_t &       filter_vectors) // <- out-param
{
    bool has_AF = false;

    /* filtered alleles R */
    filter_vectors.R.clear();
    filter_vectors.R.resize(n_alts + 1);
    filter_vectors.R[0] = 0;
    for (auto && [id, value] : record_info)
    {
        if (id == "AF")
        {
            if (!std::holds_alternative<std::vector<float>>(value))
            {
                throw decovar_error{
                  fmt::format("[Record no: {}] AF field of multi-allelic record wasn't "
                              "vector<float>.",
                              record_no)};
            }

            std::vector<float> const & afs = std::get<std::vector<float>>(value);
            if (afs.size() != n_alts)
            {
                throw decovar_error{
                  fmt::format("[Record no: {}] AF field of multi-allelic record has wrong "
                              "size: {}, but {} was expected.",
                              record_no,
                              afs.size(),
                              n_alts)};
            }

            for (size_t i = 0; i < n_alts; ++i)
                if (afs[i] < opts.rare_af_threshold)
                    filter_vectors.R[i + 1] = 1; // filtered_alleles[0] is REF and never filtered

            has_AF = true;
            break;
        }
    }

    if (!has_AF)
    {
        //TODO look for AC and AN and update those?
        throw decovar_error{fmt::format("[Record no: {}] no AF field in record.", record_no)};
    }

    /* filtered alleles A */
    filter_vectors.A.resize(n_alts);
    std::ranges::copy(filter_vectors.R.begin() + 1, filter_vectors.R.end(), filter_vectors.A.begin());

    /* filtered alleles G */
    filter_vectors.G.resize(bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1);
    for (size_t b = 0; b <= n_alts; ++b)
    {
        for (size_t a = 0; a <= b; ++a)
        {
            assert(bio::io::var::detail::vcf_gt_formula(a, b) < filter_vectors.G.size());
            filter_vectors.G[bio::io::var::detail::vcf_gt_formula(a, b)] = filter_vectors.R[a] || filter_vectors.R[b];
        }
    }
}

template <typename T>
void remove_by_indexes(std::vector<T> & vec, std::span<int const> const filter_vector)
{
    auto pred = [&](T const & elem) -> bool
    {
        ptrdiff_t i = (&elem - &*vec.begin()) % filter_vector.size(); // modulo down-maps for concat's inner vector
        assert(i >= 0);
        return filter_vector[i] != 0;
    };
    auto removed_range = std::ranges::remove_if(vec.begin(), vec.end(), pred);
    vec.resize(vec.size() - removed_range.size());
}

void update_infos(record_t::info_t &       record_info, //← in-out parameter
                  header_t const &         hdr,
                  size_t const             record_no,
                  filter_vectors_t const & filter_vectors)
{
    std::span<int const> selected_filter_vector{};

    for (auto && [_id, value] : record_info)
    {
        std::string_view         id   = _id;
        size_t const             i    = hdr.string_to_info_pos().at(id);
        header_t::info_t const & info = hdr.infos[i];

        auto visitor = bio::meta::overloaded{
          [&](auto) {
              throw decovar_error{
                fmt::format("[Record no: {}] Expected a vector when trimming field {}.", record_no, id)};
          },
          [&]<typename T>(std::vector<T> & vec)
          {
              assert(vec.size() == selected_filter_vector.size());
              remove_by_indexes(vec, selected_filter_vector);
#ifndef NDEBUG
              size_t ones = std::accumulate(selected_filter_vector.begin(), selected_filter_vector.end(), 0ul);
#endif
              assert(vec.size() == selected_filter_vector.size() - ones);
          }};

        switch (info.number)
        {
            case bio::io::var::header_number::R:
                selected_filter_vector = filter_vectors.R;
                std::visit(visitor, value);
                break;
            case bio::io::var::header_number::A:
                selected_filter_vector = filter_vectors.A;
                std::visit(visitor, value);
                break;
            default:
                break;
        }
    }
}

void update_genotypes(record_t::genotypes_t &  record_genotypes, //← in-out parameter
                      header_t const &         hdr,
                      size_t const             record_no,
                      filter_vectors_t const & filter_vectors)
{
    std::span<int const> selected_filter_vector{};

    for (auto && [_id, value] : record_genotypes)
    {
        std::string_view           id     = _id;
        size_t const               i      = hdr.string_to_format_pos().at(id);
        header_t::format_t const & format = hdr.formats[i];

        auto visitor = bio::meta::overloaded{
          [&](auto) {
              throw decovar_error{
                fmt::format("[Record no: {}] Expected a vector when trimming field {}.", record_no, id)};
          },
          [&]<typename T>(bio::ranges::concatenated_sequences<std::vector<T>> & vec)
          {
              size_t n_samples        = vec.size();
              size_t n_alleles_before = selected_filter_vector.size();
              size_t n_alleles_after =
                n_alleles_before - std::accumulate(selected_filter_vector.begin(), selected_filter_vector.end(), 0ul);

              std::pair raw_data = vec.raw_data();

              assert(vec.concat_size() == n_samples * n_alleles_before);
              assert(raw_data.second.back() == vec.concat_size());

              remove_by_indexes(raw_data.first, selected_filter_vector);
              assert(vec.concat_size() == n_samples * n_alleles_after);

              assert(vec.size() == n_samples);                 // size remains unchanged
              assert(raw_data.second.size() == n_samples + 1); // size remains unchanged

              for (size_t i = 0; i < n_samples + 1; ++i)
              {
                  assert(raw_data.second[i] == i * n_alleles_before);
                  raw_data.second[i] = i * n_alleles_after;
              }
              assert(raw_data.second.back() == raw_data.first.size());
          }};

        switch (format.number)
        {
            case bio::io::var::header_number::R:
                selected_filter_vector = filter_vectors.R;
                std::visit(visitor, value);
                break;
            case bio::io::var::header_number::A:
                selected_filter_vector = filter_vectors.A;
                std::visit(visitor, value);
                break;
            case bio::io::var::header_number::G:
                selected_filter_vector = filter_vectors.G;
                std::visit(visitor, value);
                break;
            default:
                break;
        }
    }
}

void decovar(program_options const & opts)
{
    using record_t = bio::io::var::record_default;
    bio::io::var::reader_options reader_opts{.record = record_t{}};

    bio::io::var::reader reader = (opts.input_file == "-") ? bio::io::var::reader{std::cin, bio::io::vcf{}, reader_opts}
                                                           : bio::io::var::reader{opts.input_file, reader_opts};
    bio::io::var::header const & hdr = reader.header();

    bio::io::var::writer writer = (opts.output_file == "-") ? bio::io::var::writer{std::cout, bio::io::vcf{}}
                                                            : bio::io::var::writer{opts.output_file};

    filter_vectors_t filter_vectors;
    size_t           record_no = -1;
    for (record_t & record : reader)
    {
        ++record_no;

        if (size_t n_alts = record.alt.size(); n_alts > 1ul) // multi-allelic
        {
            std::cerr << std::flush;
            std::cout << std::flush;
            fmt::print(stderr, "muli-allelic record\n");

            determine_filtered_alleles(record.info, record_no, n_alts, opts, filter_vectors);

            fmt::print(stderr, "filter_vector.A: {}\n", filter_vectors.A);
            fmt::print(stderr, "filter_vector.G: {}\n", filter_vectors.G);

            if (std::ranges::all_of(filter_vectors.A, std::identity{}))
            {
                fmt::print(stderr, "removing record {}\n", record_no);
                continue; // all alleles are removed → skip this record
            }

            // DEBUG
            fmt::print(stderr, "\n==OLD RECORD\n\n", record_no);
            std::cerr << std::flush;
            writer.push_back(record);
            // DEBUG

            /* update alts */
            remove_by_indexes(record.alt, filter_vectors.A);

            /* update info */
            update_infos(record.info, hdr, record_no, filter_vectors);

            /* update genotypes */
            update_genotypes(record.genotypes, hdr, record_no, filter_vectors);

            // DEBUG
            std::cout << std::flush;
            fmt::print(stderr, "\n==NEW RECORD\n\n", record_no);
            writer.push_back(record);
            // DEBUG
        }

        std::cerr << std::flush;
        // writer.push_back(record);
    }
}

int main(int argc, char ** argv)
{
    std::ios::sync_with_stdio(false);

#ifndef NDEBUG
    program_options o = parse_options(argc, argv);
    decovar(o);
#else
    try
    {
        program_options o = parse_options(argc, argv);
        decovar(o);
    }
    catch (sharg::parser_error const & ext)
    {
        fmt::print(stderr, "[PARSER ERROR] {}\n", ext.what());
        return -1;
    }
    catch (bio::io::bio_error const & ext)
    {
        fmt::print(stderr, "[BioC++ I/O error] {}\n", ext.what());
        return -1;
    }
    catch (decovar_error const & ext)
    {
        fmt::print(stderr, "[deCoVar error] {}\n", ext.what());
        return -1;
    }
#endif
    return 0;
}
