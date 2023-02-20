// MIT License
//
// Copyright (c) 2023 deCODE Genetics
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "binalleles.hpp"

#include <bits/ranges_algo.h>
#include <cstddef>
#include <cstdint>
#include <strings.h>
#include <variant>

#include <bio/io/exception.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/misc.hpp>
#include <bio/io/var/reader.hpp>
#include <bio/io/var/writer.hpp>
#include <bio/ranges/to.hpp>
#include <bio/ranges/views/zip.hpp>

#include <sharg/all.hpp>

#include "../generator.hpp"
#include "../misc.hpp"
#include "bio/io/misc.hpp"
#include "bio/io/var/record.hpp"

namespace _binalleles
{

program_options parse_options(sharg::parser & parser)
{
    program_options opts;

    parser.add_flag(
      opts.verbose,
      sharg::config{.short_id = 'v', .long_id = "verbose", .description = "Print diagnostics to stderr."});

    parser.add_subsection("Input / Output:");
    parser.add_positional_option(opts.input_file,
                                 sharg::config{.description = "Path to input file or '-' for stdin.",
                                               .required    = true,
                                               .validator   = input_file_or_stdin_validator{{"vcf", "vcf.gz", "bcf"}}});
    parser.add_option(opts.output_file,
                      sharg::config{
                        .short_id    = 'o',
                        .long_id     = "output",
                        .description = "Path to output file or '-' for stdout.",
                        .validator   = output_file_or_stdout_validator{sharg::output_file_open_options::create_new,
                                                                       {"vcf", "vcf.gz", "bcf"}}
    });

    parser.add_option(opts.output_file_type,
                      sharg::config{
                        .short_id    = 'O',
                        .long_id     = "output-type",
                        .description = "Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), "
                                       "uncompressed VCF (v); or use automatic (a) detection. Use the -Ou option "
                                       "when piping between subcommands to speed up performance by removing "
                                       "unnecessary compression/decompression and VCF←→BCF conversion.",
                        .validator   = sharg::value_list_validator{'a', 'b', 'u', 'z', 'v'}
    });

    parser.add_subsection("Allele binning by length:");

    parser.add_line(
      "Splits every n-allelic record into up to (n-1) records with two pseudo-alleles each. The (pseudo) REF allele "
      "encompasses all original alleles shorter than the threshold; and the (pseudo) ALT allele encompasses "
      "all original alleles whose length is >= the threshold. The threshold changes within each batch, so that "
      "the first record has only the shortest allele in REF and all others in ALT. Subsequently, the next shortest "
      "allele is moved from the ALT bin into the REF bin; until the ALT bin only contains the longest allele.",
      true);

    parser.add_flag(opts.bin_by_length,
                      sharg::config{
                        .long_id     = "bin-by-length",
                        .description = "Activates this option."
    });

    parser.add_flag(opts.same_length_splits,
                      sharg::config{
                        .long_id     = "same-length-splits",
                        .description = "By default, records are skipped where the split happens between alleles of "
                        "the same length. This options enables writing of all records."
    });

    parser.add_subsection("Performance:");
    parser.add_option(opts.threads,
                      sharg::config{
                        .short_id    = '@',
                        .long_id     = "threads",
                        .description = "Maximum number of threads to use.",
                        .validator   = sharg::arithmetic_range_validator{2u, std::thread::hardware_concurrency() * 2}
    });

    parser.parse();
    return opts;
}

template <typename int_t>
bio::ranges::concatenated_sequences<std::vector<int_t>> & establish_PLs(
    bio::io::var::genotype_element_value_type<bio::io::ownership::deep> & variant,
    size_t const n_samples)
{
    if (!std::holds_alternative<bio::ranges::concatenated_sequences<std::vector<int_t>>>(variant))
    {
        bio::ranges::concatenated_sequences<std::vector<int_t>> newGTs;
        concatenated_sequences_create_scaffold(newGTs, n_samples, 3ul);
        variant = std::move(newGTs);
    }
    return std::get<bio::ranges::concatenated_sequences<std::vector<int_t>>>(variant);
}

void main(sharg::parser & parser)
{
    program_options opts = parse_options(parser);

    size_t threads        = opts.threads - 1; // subtract one for the main thread
    size_t reader_threads = threads / 3;
    size_t writer_threads = threads - reader_threads;

    /* setup reader */
    bio::io::var::reader_options reader_opts{.record = record_t{},
                                             .stream_options =
                                               bio::io::transparent_istream_options{.threads = reader_threads + 1}};
    bio::io::var::reader reader = (opts.input_file == "-") ? bio::io::var::reader{std::cin, bio::io::vcf{}, reader_opts}
                                                           : bio::io::var::reader{opts.input_file, reader_opts};

    /* setup writer */
    bio::io::var::writer writer = create_writer(opts.output_file, opts.output_file_type, writer_threads);

    /* ========= setup header =========== */
    if (opts.bin_by_length) // we need to create a new header
    {
        bio::io::var::header new_hdr = reader.header(); // create copy

        new_hdr.infos.clear();
        bio::io::var::header::format_t ref{
            .id = "REFBIN_INDEXES",
            .number = bio::io::var::header_number::dot,
            .type = "Integer",
            .type_id = bio::io::var::value_type_id::vector_of_int16, //TODO configure size?
            .description = "Indexes of original alleles binned as the reference.",
        };
        new_hdr.infos.push_back(std::move(ref));

        bio::io::var::header::format_t ref_max{
            .id = "REFBIN_MAXLEN",
            .number = 1,
            .type = "Integer",
            .type_id = bio::io::var::value_type_id::int32,
            .description = "Maximum allele length in REFBIN.",
        };
        new_hdr.infos.push_back(std::move(ref_max));


        bio::io::var::header::format_t alt{
            .id = "ALTBIN_INDEXES",
            .number = bio::io::var::header_number::dot,
            .type = "Integer",
            .type_id = bio::io::var::value_type_id::vector_of_int16, //TODO configure size?
            .description = "Indexes of original alleles binned as the ALT.",
        };
        new_hdr.infos.push_back(std::move(alt));

        bio::io::var::header::format_t alt_min{
            .id = "ALTBIN_MINLEN",
            .number = 1,
            .type = "Integer",
            .type_id = bio::io::var::value_type_id::int32,
            .description = "Minimum allale length in ALTBIN.",
        };
        new_hdr.infos.push_back(std::move(alt_min));


        new_hdr.formats.clear();
        new_hdr.formats.push_back(bio::io::var::reserved_formats.at("GT"));
        new_hdr.formats.push_back(bio::io::var::reserved_formats.at("PL"));

        new_hdr.add_missing();
        writer.set_header(std::move(new_hdr));
    }
    else // just use existing header as-is
    {
        writer.set_header(reader.header());
    }

    bio::io::var::header const & hdr = writer.header();
    if (hdr.column_labels.size() < 10)
        throw decovar_error{"VCF file contains no samples."};
    size_t const n_samples = hdr.column_labels.size() - 9;

    /* caches */
    size_t   record_no = -1; // always refers to #record in input even if more records are created
    /* temporary record – do we need dynamic int-width?*/
    record_t new_rec;
    new_rec.alt = {".", "."};
    new_rec.info.emplace_back("REFBIN_MAXLEN", int32_t{});
    new_rec.info.emplace_back("ALTBIN_MINLEN", int32_t{});
    new_rec.info.emplace_back("REFBIN_INDEXES", std::vector<int16_t>{});
    new_rec.info.emplace_back("ALTBIN_INDEXES", std::vector<int16_t>{});
    new_rec.genotypes.emplace_back("GT", std::vector<std::string>{});
    using pl_t = bio::ranges::concatenated_sequences<std::vector<int16_t>>;
    new_rec.genotypes.emplace_back("PL", pl_t{});

    /* shortcuts into the variants */
    int32_t &                  refbin_max     = std::get<int32_t>(new_rec.info[0].value);
    int32_t &                  altbin_min     = std::get<int32_t>(new_rec.info[1].value);
    std::vector<int16_t> &     refbin_indexes = std::get<std::vector<int16_t>>(new_rec.info[2].value);
    std::vector<int16_t> &     altbin_indexes = std::get<std::vector<int16_t>>(new_rec.info[3].value);
    std::vector<std::string> & out_GTs        = std::get<std::vector<std::string>>(new_rec.genotypes[0].value);
    out_GTs.resize(n_samples);

    /* ========= define steps =========== */

    /* pre */
    auto pre_fn = [&](record_t & record) -> record_t &
    {
        ++record_no;
        return record;
    };
    auto pre_view = std::views::transform(pre_fn);

    /* remove rare alleles */
    auto bin_by_length_fn = [&](record_t & record) -> std::generator<record_t &>
    {
        size_t const n_alts    = record.alt.size();
        size_t const n_alleles = n_alts + 1;

        if (n_alts <= 1ul || !opts.bin_by_length /* || !record.genotypes.contains("GT")*/)
        {
            co_yield record;
            co_return;
        }

        /* this crap can go once we have dictionaries */
        {
            bool has_PL = false;
            for (auto && [key, value] : record.genotypes)
                    if (key == "PL")
                        has_PL = true;
            if (!has_PL)
            {
                co_yield record;
                co_return;
            }
        }

        //                   length, index
        std::vector<std::pair<size_t, size_t>> allele_lengths;
        allele_lengths.resize(n_alleles);
        allele_lengths[0] = {record.ref.size(), 0};

        auto lengths_v = allele_lengths | std::views::elements<0>;
        std::ranges::copy(record.alt | std::views::transform(std::ranges::size), lengths_v.begin() + 1);

        auto indexes_v = allele_lengths | std::views::elements<1>;
        std::ranges::copy(std::views::iota(1ul, n_alleles), indexes_v.begin() + 1);

        // fmt::print(stderr, "allele_lengths: {}\n", allele_lengths);
        std::ranges::sort(allele_lengths);
        // fmt::print(stderr, "allele_lengths (sorted): {}\n", allele_lengths);

        new_rec.chrom = record.chrom;
        new_rec.pos  = record.pos;


        for (size_t i = 0; i < indexes_v.size() - 1; ++i)
        {
            refbin_max = lengths_v[i];
            altbin_min = lengths_v[i + 1];

            if (refbin_max == altbin_min && !opts.same_length_splits) // lengths shall not be present in both groups
                continue;

            if (record.id != ".")
                new_rec.id = fmt::format("{}_div_{}", record.id, i);

            refbin_indexes.clear();
            std::ranges::copy(indexes_v | std::views::take(i + 1), std::back_inserter(refbin_indexes));

            altbin_indexes.clear();
            std::ranges::copy(indexes_v | std::views::drop(i + 1), std::back_inserter(altbin_indexes));

            auto visitor = bio::meta::overloaded{
                [] (auto const &) { throw decovar_error{"PL field was in wrong state"}; },
                [&] <typename int_t> (bio::ranges::concatenated_sequences<std::vector<int_t>> const & in_PLs)
            {
                assert(n_samples == in_PLs.size());
                if (in_PLs.concat_size() != n_samples * (bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1))
                {
                    throw decovar_error{
                    "[Record no: {}] Currently, every sample must be diploid and must contain the "
                    "full number of PL values (e.g. no single '.' placeholder allowed).",
                    record_no};
                }

                bio::ranges::concatenated_sequences<std::vector<int_t>> & out_PLs =
                    establish_PLs<int_t>(new_rec.genotypes[1].value, n_samples);

                for (size_t j = 0; j < n_samples; ++j)
                {
                    std::span<int_t const> const in_PL  = in_PLs[j];
                    assert(in_PL.size() == bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1);
                    std::span<int_t> const       out_PL = out_PLs[j];

                    /* 0/0 value */
                    out_PL[0] = std::numeric_limits<int_t>::max();
                    for (size_t const b : refbin_indexes)
                        for (size_t const a : refbin_indexes)
                            if (a <= b)
                                out_PL[0] = std::min<int_t>(out_PL[0], in_PL[bio::io::var::detail::vcf_gt_formula(a, b)]);

                    /* 0/1 value */
                    out_PL[1] = std::numeric_limits<int_t>::max();
                    for (size_t const b : refbin_indexes)
                        for (size_t const a : altbin_indexes)
                            if (a <= b)
                                out_PL[1] = std::min<int_t>(out_PL[1], in_PL[bio::io::var::detail::vcf_gt_formula(a, b)]);
                    for (size_t const b : altbin_indexes)
                        for (size_t const a : refbin_indexes)
                            if (a <= b)
                                out_PL[1] = std::min<int_t>(out_PL[1], in_PL[bio::io::var::detail::vcf_gt_formula(a, b)]);

                    /* 1/1 value */
                    out_PL[2] = std::numeric_limits<int_t>::max();
                    for (size_t const b : altbin_indexes)
                        for (size_t const a : altbin_indexes)
                            if (a <= b)
                                out_PL[2] = std::min<int_t>(out_PL[2], in_PL[bio::io::var::detail::vcf_gt_formula(a, b)]);

                    switch(std::ranges::min_element(out_PL) - out_PL.begin())
                    {
                        case 0:
                            out_GTs[j] = "0/0";
                            break;
                        case 1:
                            out_GTs[j] = "0/1";
                            break;
                        case 2:
                            out_GTs[j] = "1/1";
                            break;
                        default:
                            BIOCPP_UNREACHABLE;
                            break;
                    }
                }
            }};

            for (auto && [key, value] : record.genotypes)
                if (key == "PL")
                    std::visit(visitor, value), ({break;});

            co_yield new_rec;

        }
    };
    auto bin_by_length_view = std::views::transform(bin_by_length_fn) | views_cojoin;

    /* ========= create and execute pipeline =========== */
    reader | pre_view | bin_by_length_view | writer;
}

}
