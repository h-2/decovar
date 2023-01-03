#include "bio/io/exception.hpp"
#include "bio/io/var/header.hpp"
#include "sharg/validators.hpp"
#include <bio/io/var/reader.hpp>
#include <bio/io/var/writer.hpp>
#include <cstddef>
#include <sharg/all.hpp>
#include <variant>

using record_t = bio::io::var::record_default;
using header_t = bio::io::var::header;

struct decovar_error : std::runtime_error
{
    decovar_error(std::string msg) : std::runtime_error{msg} {}
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
    parser.info.author = "Hannes Hauswedell";
    parser.info.short_description = "Reduce allele complexity in VCF files.";
    parser.info.version = "0.1.0";
    parser.add_positional_option(o.input_file,
                                 sharg::config{.description = "Path to input file or '-' for stdin.",
                                               .required = true,
                                               .validator = sharg::input_file_validator{{"vcf", "vcf.gz", "bcf"}}});
    parser.add_option(o.output_file,
        sharg::config{.short_id = 'o',
                      .long_id = "output",
                      .description = "Path to output file or '-' for stdout.",
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::create_new,
                                                                {"vcf", "vcf.gz", "bcf"}}});

    parser.add_option(o.rare_af_threshold,
    sharg::config{.long_id = "rare-af-thresh",
                  .description = "Threshold below which alleles are considered rare.",
                  .validator = sharg::arithmetic_range_validator{1e-10, 1e-1}});

    parser.parse();
    return o;
}


void determine_filtered_alleles(record_t::info_t const & record_info,
                                size_t const record_no,
                                size_t const R,
                                program_options const & opts,
                                std::vector<int> & filtered_alleles) // <- out-param
{
    bool has_AF = false;

    for (auto && [ id, value ] : record_info)
    {
        if (id == "AF")
        {
            if (!std::holds_alternative<std::vector<float>>(value))
            {
                throw decovar_error{fmt::format("[Record no: {}] AF field of multi-allelic record wasn't "
                                                "vector<float>.", record_no)};
            }

            std::vector<float> const & afs = std::get<std::vector<float>>(value);
            if (afs.size() != R)
            {
                throw decovar_error{fmt::format("[Record no: {}] AF field of multi-allelic record has wrong "
                                                "size: {}, but {} was expected.", record_no, afs.size(), R)};
            }

            for (size_t i = 0; i < R; ++i)
                if (afs[i] < opts.rare_af_threshold)
                    filtered_alleles[i] = 1;

            has_AF = true;
            break;
        }
    }

    if (!has_AF)
    {
        //TODO look for AC and AN and update those?
        throw decovar_error{fmt::format("[Record no: {}] no AF field in record.", record_no)};
    }
}

template <typename T>
void remove_by_indexes(std::vector<T> & vec,
                       bool const skip_vector_1st_elem,
                       std::vector<int> const & filter_vector)
{
    assert(vec.size() - skip_vector_1st_elem == filter_vector.size());

    auto const beg = vec.begin() + skip_vector_1st_elem;

    auto pred = [&] (T const & elem) -> bool
    {
        ptrdiff_t i = &elem - &*beg;
        assert(i >= 0);
        return filter_vector[i];
    };
    auto ret_range = std::ranges::remove_if(beg, vec.end(), pred);
    size_t new_size = ret_range.size() + skip_vector_1st_elem;

#ifndef NDEBUG
    size_t ones = std::accumulate(filter_vector.begin(), filter_vector.end(), 0ul);
    assert(ones + skip_vector_1st_elem == new_size);
#endif

    vec.resize(new_size);
}

void update_infos(record_t::info_t & record_info, //← in-out parameter
                  header_t const & hdr,
                  size_t const record_no,
                  std::vector<int> const filtered_alleles)
{
    for (auto && [ _id, value ] : record_info)
    {
        std::string_view id = _id;
        size_t const i = hdr.string_to_info_pos().at(id);
        header_t::info_t const & info = hdr.infos[i];

        bool is_R = false;
        auto visitor = bio::meta::overloaded{
            [&] (auto)
            {
                throw decovar_error{fmt::format("[Record no: {}] Expected a vector when trimming field {}.",
                                                record_no,
                                                id)};
            },
            [&] <typename T> (std::vector<T> & vec)
            {
                remove_by_indexes(vec, is_R, filtered_alleles);
            }
        };

        switch(info.number)
        {
            case bio::io::var::header_number::R:
                is_R = true;
                [[fallthrough]];
            case bio::io::var::header_number::A:
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

    bio::io::var::reader reader = (opts.input_file == "-") ? bio::io::var::reader{std::cin, bio::io::vcf{}, reader_opts} :
    bio::io::var::reader{opts.input_file, reader_opts};
    bio::io::var::header const & hdr = reader.header();

    bio::io::var::writer writer = (opts.output_file == "-") ? bio::io::var::writer{std::cout, bio::io::vcf{}} : bio::io::var::writer{opts.output_file};

    std::vector<int/*bool*/> filtered_alleles; // bool vector of allele number with filtered yes/no
    size_t record_no = -1;
    for (record_t & record : reader)
    {
        ++record_no;

        record.genotypes.clear(); // improve readability temporarily

        if (size_t R = record.alt.size(); R > 1ul) // multi-allelic
        {
            std::cerr << std::flush;
            std::cout << std::flush;
            fmt::print(stderr, "muli-allelic record\n");

            filtered_alleles.clear();
            filtered_alleles.resize(R);
            determine_filtered_alleles(record.info, record_no, R, opts, filtered_alleles);

            fmt::print(stderr, "filtered_alleles: {}\n", filtered_alleles);

            if (std::ranges::all_of(filtered_alleles, std::identity{}))
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
            remove_by_indexes(record.alt, false, filtered_alleles);

            /* update info */
            update_infos(record.info, hdr, record_no, filtered_alleles);

            /* update genotypes */
            // TODO

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



int main(int argc, char **argv)
{
    std::ios::sync_with_stdio(false);

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

    return 0;
}
