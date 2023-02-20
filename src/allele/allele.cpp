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

#include "../generator.hpp"
#include "../misc.hpp"
#include "localise.hpp"
#include "remove.hpp"
#include "split.hpp"

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
    parser.add_line(
      "Allows removing certain alleles from multi-allelic records. All fields with A, R or G multiplicity"
      " have the respective elements removed. The GT field is updated to contain the new indexes.",
      true);

    parser.add_option(opts.rare_af_threshold,
                      sharg::config{
                        .long_id     = "rare-af-thresh",
                        .description = "For multi-allelic records, remove alleles with AF < than threshold. "
                                       "0 → remove none.",
                        .validator   = sharg::arithmetic_range_validator{0.0, 1.0}
    });

    parser.add_subsection("Divide alleles into multiple records:");

    parser.add_line("Multi-allelic records are split into two or more with some alleles each.", true);

    parser.add_option(opts.split_by_length,
                      sharg::config{
                        .long_id     = "split-by-length",
                        .description = "Alleles shorter than this will stay in this record; longer "
                                       "ones are moved into a separate one. 0 → no splitting.",
                        .validator   = sharg::arithmetic_range_validator{0ul, 100'000ul}
    });

    parser.add_subsection("Allele localisation:");

    parser.add_line(
      "Determine the \"locally relevant\" alleles per sample (using the PL-field) and store their "
      "indexes in the newly added LAA-field. Note that the reference allele "
      "denoted by 0 is always considered locally relevant without being listed in LAA.\n"
      "The PL field and AD field are then renamed to LPL and LAD and subsampled to only contain "
      "information for the local alleles.",
      true);

    parser.add_option(opts.local_alleles,
                      sharg::config{
                        .short_id    = 'L',
                        .long_id     = "local-alleles",
                        .description = "For multi-allelic records with more than L alleles, transform global alleles "
                                       "to local alleles. 0 → never transform.",
                        .validator   = sharg::arithmetic_range_validator{0, 127}
    });

    parser.add_option(opts.keep_global_fields,
                      sharg::config{.long_id     = "keep-global-fields",
                                    .description = "If set, PL and AD fields are kept in addition to LPL and LAD."});

    parser.add_option(opts.transform_all,
                      sharg::config{.long_id     = "transform-all",
                                    .description = "If set, records with fewer than L alleles will still get an "
                                                   "LAA-field and have their PL/AD renamed to LPL/LAD. This increases "
                                                   "file size and provides no advantage other than enabling same "
                                                   "FORMATs for all records."});

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
    bio::io::var::writer writer = create_writer(opts.output_file, opts.output_file_type, writer_threads);

    /* ========= setup header =========== */
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
    size_t             record_no = -1; // always refers to #record in input even if more records are created
    _remove::cache_t   filter_vectors;
    _localise::cache_t localise_cache;

    /* ========= define steps =========== */

    /* pre */
    auto pre_fn = [&](record_t & record) -> record_t &
    {
        ++record_no;
        return record;
    };
    auto pre_view = std::views::transform(pre_fn);

    /* remove rare alleles */
    auto remove_rare_alleles_fn = [&](record_t & record) -> std::generator<record_t &>
    {
        if (record.alt.size() > 1ul && opts.rare_af_threshold != 0.0)
        {
            log(opts, "↓ record no {} allelle-removal begin.\n", record_no);
            bool all_alleles_removed = _remove::remove_rare_alleles(record, record_no, hdr, opts, filter_vectors);
            log(opts, "↑ record no {} allelle-removal end.\n", record_no);
            if (all_alleles_removed)
                co_return;
        }

        co_yield record;
    };
    auto remove_rare_alleles_view = std::views::transform(remove_rare_alleles_fn) | views_cojoin;

    /* split */
    auto split_fn = [&](record_t & record) -> std::generator<record_t &>
    {
        if (opts.split_by_length > 0 && _split::needs_splitting(record, opts))
        {
            log(opts, "↓ record no {} splitting-by-length begin.\n", record_no);

            /* create second record */
            record_t record0 = record;
            if (record0.id != ".")
            {
                record0.id += "_split1";
                record.id += "_split2";
            }

            // short alleles (remove gt)
            _split::remove_alleles(record0, record_no, _split::gt, hdr, opts, filter_vectors);

            // long alleles (remove leq)
            _split::remove_alleles(record, record_no, _split::leq, hdr, opts, filter_vectors);

            log(opts, "↑ record no {} splitting-by-length end.\n", record_no);

            co_yield record0;
        }

        co_yield record;
    };
    auto split_view = std::views::transform(split_fn) | views_cojoin;

    /* localise */
    auto localise_fn = [&](record_t & record) -> record_t &
    {
        if (opts.local_alleles != 0)
        {
            if (record.alt.size() > opts.local_alleles)
            {
                log(opts, "↓ record no {} allelle-localisation begin.\n", record_no);
                _localise::localise_alleles(record, record_no, hdr, opts, localise_cache);
                log(opts, "↑ record no {} allelle-localisation end.\n", record_no);
            }
            else if (opts.transform_all)
            {
                log(opts, "↓ record no {} allelle-pseudo-localisation begin.\n", record_no);
                _localise::pseudo_localise_alleles(record, record_no, hdr, opts, localise_cache);
                log(opts, "↑ record no {} allelle-pseudo-localisation end.\n", record_no);
            }
        }

        return record;
    };
    auto localise_view = std::views::transform(localise_fn);

    /* ========= create pipeline =========== */
    auto pipeline = reader | pre_view | remove_rare_alleles_view | split_view | localise_view;

    /* ========= iterate =========== */
    for (record_t & record : pipeline)
    {
        /* finally write the (modified) record */
        writer.push_back(record);

        /* salvage memory */
        if (opts.local_alleles != 0 && ((record.alt.size() > opts.local_alleles) || opts.transform_all))
            _localise::salvage_cache(record, localise_cache);
    }
}
