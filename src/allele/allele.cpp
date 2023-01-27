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

#include "allele.hpp"

#include <cstddef>
#include <variant>

#include <bio/io/exception.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/misc.hpp>
#include <bio/io/var/reader.hpp>
#include <bio/io/var/writer.hpp>

#include <sharg/all.hpp>

#include "../misc.hpp"
#include "localise.hpp"
#include "remove_rare.hpp"

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

    parser.add_subsection("Algorithm:");
    parser.add_option(opts.rare_af_threshold,
                      sharg::config{
                        .long_id     = "rare-af-thresh",
                        .description = "For multi-allelic records, remove alleles with AF < than threshold. "
                                       "0 → remove none.",
                        .validator   = sharg::arithmetic_range_validator{0.0, 1.0}
    });

    parser.add_option(opts.local_alleles,
                      sharg::config{
                        .short_id    = 'L',
                        .long_id     = "local-alleles",
                        .description = "For multi-allelic records with more than L alleles, transform global alleles "
                                       "to local alleles. 0 → never transform.",
                        .validator   = sharg::arithmetic_range_validator{0, 127}
    });

    parser.add_option(opts.remove_global_alleles,
                      sharg::config{.long_id     = "remove-global-alleles",
                                    .description = "Remove global alleles for whome local alleles have been created."});

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

void allele(sharg::parser & parser)
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
    bio::io::var::writer_options writer_opts{.stream_options =
                                               bio::io::transparent_ostream_options{.threads = writer_threads + 1}};

    bio::io::var::writer writer = (opts.output_file == "-")
                                    ? bio::io::var::writer{std::cout, bio::io::vcf{}, writer_opts}
                                    : bio::io::var::writer{opts.output_file, writer_opts};

    /* setup header */
    if (opts.local_alleles > 0ul) // we need to create a new header
    {
        bio::io::var::header new_hdr = reader.header(); // create copy

        if (!new_hdr.string_to_format_pos().contains("LAA"))
            new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LAA"));
        if (new_hdr.string_to_format_pos().contains("AD") && !new_hdr.string_to_format_pos().contains("LAD"))
            new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LAD"));
        // if (new_hdr.string_to_format_pos().contains("GT") && !new_hdr.string_to_format_pos().contains("LGT"))
        //     new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LGT"));
        if (new_hdr.string_to_format_pos().contains("PL") && !new_hdr.string_to_format_pos().contains("LPL"))
            new_hdr.formats.push_back(bio::io::var::reserved_formats.at("LPL"));

        new_hdr.add_missing();
        writer.set_header(std::move(new_hdr));
    }
    else // just use existing header as-is
    {
        writer.set_header(reader.header());
    }

    bio::io::var::header const & hdr = writer.header();

    /* caches */
    size_t           record_no = -1;
    filter_vectors_t filter_vectors;
    localise_cache_t localise_cache;

    /* iterate */
    for (record_t & record : reader)
    {
        ++record_no;

        /* remove rare alleles */
        if (record.alt.size() > 1ul && opts.rare_af_threshold != 0.0)
        {
            log(opts, "↓ record no {} allelle-removal begin.\n", record_no);
            bool all_alleles_removed = remove_rare_alleles(record, record_no, hdr, opts, filter_vectors);
            log(opts, "↑ record no {} allelle-removal end.\n", record_no);
            if (all_alleles_removed == true)
                continue;
        }

        /* localise alleles */
        if (record.alt.size() > opts.local_alleles && opts.local_alleles != 0)
        {
            log(opts, "↓ record no {} allelle-localisation begin.\n", record_no);
            localise_alleles(record, record_no, hdr, opts, localise_cache);
            log(opts, "↑ record no {} allelle-localisation end.\n", record_no);
        }

        /* finally write the (modified) record */
        writer.push_back(record);

        /* salvage memory */
        if (record.alt.size() > opts.local_alleles && opts.local_alleles != 0)
            salvage_localise_cache(record, localise_cache);
    }
}
