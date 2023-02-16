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

    parser.add_subsection("Remove rare alleles:");

    parser.add_flag(opts.bin_by_length,
                      sharg::config{
                        .long_id     = "bin-by-length",
                        .description = "TODO"
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
    // if (opts.bin_by_length > 0ul) // we need to create a new header
    // {
    //     bio::io::var::header new_hdr = reader.header(); // create copy
    //
    //     if (!new_hdr.string_to_format_pos().contains("LAA"))
    //         new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LAA"));
    //     if (new_hdr.string_to_format_pos().contains("AD") && !new_hdr.string_to_format_pos().contains("LAD"))
    //         new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LAD"));
    //     // if (new_hdr.string_to_format_pos().contains("GT") && !new_hdr.string_to_format_pos().contains("LGT"))
    //     //     new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LGT"));
    //     if (new_hdr.string_to_format_pos().contains("PL") && !new_hdr.string_to_format_pos().contains("LPL"))
    //         new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LPL"));
    //
    //     new_hdr.add_missing();
    //     writer.set_header(std::move(new_hdr));
    // }
    // else // just use existing header as-is
    {
        writer.set_header(reader.header());
    }

    // bio::io::var::header const & hdr = writer.header();

    /* caches */
    size_t             record_no = -1; // always refers to #record in input even if more records are created

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

        fmt::print(stderr, "allele_lengths: {}\n", allele_lengths);
        std::ranges::sort(allele_lengths);
        fmt::print(stderr, "allele_lengths (sorted): {}\n", allele_lengths);

        record_t new_rec{.chrom = record.chrom, .pos = record.pos, .alt{".", "."}};

        for (size_t i = 0; i < indexes_v.size() - 1; ++i)
        {
            if (record.id != ".")
                new_rec.id = fmt::format("{}_bin_{}", record.id, i);

            new_rec.genotypes.clear(); //TODO salvage memory

            auto group1 = indexes_v | std::views::take(i + 1);
            auto group2 = indexes_v | std::views::drop(i + 1);

            fmt::print(stderr, "group1: {}\n", group1);
            fmt::print(stderr, "group2: {}\n", group2);

            auto visitor = bio::meta::overloaded{
                [] (auto const &) { throw decovar_error{"PL field was in wrong state"}; },
                [&] <typename int_t> (bio::ranges::concatenated_sequences<std::vector<int_t>> const & in_PLs)
            {
                size_t const n_samples = in_PLs.size();
                if (in_PLs.concat_size() != n_samples * (bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1))
                {
                    throw decovar_error{
                    "[Record no: {}] Currently, every sample must be diploid and must contain the "
                    "full number of PL values (e.g. no single '.' placeholder allowed).",
                    record_no};
                }

                bio::ranges::concatenated_sequences<std::vector<int_t>> out_PLs;
                concatenated_sequences_create_scaffold(out_PLs, n_samples, 3ul);

                std::vector<std::string> out_GTs;
                out_GTs.resize(n_samples);

                for (size_t j = 0; j < n_samples; ++j)
                {
                    std::span<int_t const> const in_PL  = in_PLs[j];
                    assert(in_PL.size() == bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1);
                    std::span<int_t> const       out_PL = out_PLs[j];

                    /* 0/0 value */
                    out_PL[0] = std::numeric_limits<int_t>::max();
                    for (size_t const b : group1)
                        for (size_t const a : group1)
                            if (a <= b)
                                out_PL[0] = std::min<int_t>(out_PL[0], in_PL[bio::io::var::detail::vcf_gt_formula(a, b)]);

                    /* 0/1 value */
                    out_PL[1] = std::numeric_limits<int_t>::max();
                    for (size_t const b : group1)
                        for (size_t const a : group2)
                            if (a <= b)
                                out_PL[1] = std::min<int_t>(out_PL[1], in_PL[bio::io::var::detail::vcf_gt_formula(a, b)]);
                    for (size_t const b : group2)
                        for (size_t const a : group1)
                            if (a <= b)
                                out_PL[1] = std::min<int_t>(out_PL[1], in_PL[bio::io::var::detail::vcf_gt_formula(a, b)]);

                    /* 1/1 value */
                    out_PL[2] = std::numeric_limits<int_t>::max();
                    for (size_t const b : group2)
                        for (size_t const a : group2)
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

                new_rec.genotypes.emplace_back("GT", std::move(out_GTs));
                new_rec.genotypes.emplace_back("PL", std::move(out_PLs));
            }};

            for (auto && [key, value] : record.genotypes)
                if (key == "PL")
                    std::visit(visitor, value), ({break;});

            co_yield new_rec;

        }
    };
    auto bin_by_length_view = std::views::transform(bin_by_length_fn) | std::views::join;

    /* ========= create and execute pipeline =========== */
    reader | pre_view | bin_by_length_view | writer;
}

}
