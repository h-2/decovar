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

#include <cstddef>
#include <sharg/all.hpp>
#include <variant>

#include <bio/io/exception.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/reader.hpp>
#include <bio/io/var/writer.hpp>

#include "misc.hpp"
#include "remove_rare_alleles.hpp"

program_options parse_options(int const argc, char const * const * const argv)
{
    program_options o;

    sharg::parser parser{"deCoVar", argc, argv};
    parser.info.author            = "Hannes Hauswedell";
    parser.info.short_description = "deCODE variant tools.";
    parser.info.version           = "0.1.0";

    parser.add_flag(o.verbose,
                      sharg::config{
                        .short_id    = 'v',
                        .long_id     = "verbose",
                        .description = "Print diagnostics to stderr."
    });


    parser.add_subsection("Input / Output :");
    parser.add_positional_option(o.input_file,
                                 sharg::config{.description = "Path to input file or '-' for stdin.",
                                               .required    = true,
                                               .validator   = input_file_or_stdin_validator{{"vcf", "vcf.gz", "bcf"}}});
    parser.add_option(o.output_file,
                      sharg::config{
                        .short_id    = 'o',
                        .long_id     = "output",
                        .description = "Path to output file or '-' for stdout.",
                        .validator   = output_file_or_stdout_validator{sharg::output_file_open_options::create_new,
                                                                       {"vcf", "vcf.gz", "bcf"}}
    });

    parser.add_subsection("Algorithm:");
    parser.add_option(o.rare_af_threshold,
                      sharg::config{
                        .long_id     = "rare-af-thresh",
                        .description = "For multia-allelic records, remove alleles with AF < than threshold. "
                                       "0 → remove none.",
                        .validator   = sharg::arithmetic_range_validator{0.0, 1.0}
    });

    parser.add_subsection("Performance:");
    parser.add_option(o.threads,
                      sharg::config{
                        .short_id    = '@',
                        .long_id     = "threads",
                        .description = "Maximum number of threads to use.",
                        .validator   = sharg::arithmetic_range_validator{2u, std::thread::hardware_concurrency() * 2}
    });

    parser.parse();
    return o;
}

void decovar(program_options const & opts)
{

    size_t threads        = opts.threads - 1; // subtract one for the main thread
    size_t reader_threads = threads / 3;
    size_t writer_threads = threads - reader_threads;

    /* setup reader */
    bio::io::var::reader_options reader_opts{
        .record = record_t{},
        .stream_options = bio::io::transparent_istream_options{.threads = reader_threads + 1}
    };
    bio::io::var::reader reader = (opts.input_file == "-") ? bio::io::var::reader{std::cin, bio::io::vcf{}, reader_opts}
                                                           : bio::io::var::reader{opts.input_file, reader_opts};
    bio::io::var::header const & hdr = reader.header();

    /* setup reader */
    bio::io::var::writer_options writer_opts{
        .stream_options = bio::io::transparent_ostream_options{.threads = writer_threads + 1}
    };

    bio::io::var::writer writer = (opts.output_file == "-") ? bio::io::var::writer{std::cout, bio::io::vcf{}, writer_opts}
                                                            : bio::io::var::writer{opts.output_file, writer_opts};
    writer.set_header(hdr);

    /* caches */
    filter_vectors_t filter_vectors;
    size_t           record_no = -1;

    /* iterate */
    for (record_t & record : reader)
    {
        ++record_no;

        /* remove rare alleles */
        if (record.alt.size() > 1 && opts.rare_af_threshold > 0.0)
        {
            bool all_alleles_removed = remove_rare_alleles(record, record_no, hdr, opts, filter_vectors);
            if (all_alleles_removed == true)
                continue;
        }

        /* finally write the (modified) record */
        writer.push_back(record);
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
