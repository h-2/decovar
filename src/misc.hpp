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

#pragma once

#include <bio/alphabet/fmt.hpp>
#include <bio/io/format/vcf.hpp>
#include <bio/io/stream/compression.hpp>
#include <bio/io/var/header.hpp>
#include <bio/io/var/record.hpp>
#include <bio/io/var/writer.hpp>

#include <sharg/all.hpp>

using record_t = bio::io::var::record_default;
using header_t = bio::io::var::header;

struct decovar_error : std::runtime_error
{
    template <typename... T>
    decovar_error(fmt::format_string<T...> format, T &&... args) :
      std::runtime_error{fmt::format(format, std::forward<T>(args)...)}
    {}
};

template <typename opts_t, typename... T>
inline void log(opts_t const & opts, fmt::format_string<T...> format, T &&... args)
{
    if (opts.verbose)
    {
        fmt::print(stderr, "[decovar log] ");
        fmt::print(stderr, format, std::forward<T>(args)...);
    }
}

// ============================================================================
// Argument parsing
// ============================================================================

class input_file_or_stdin_validator : public sharg::input_file_validator
{
private:
    using base_t = sharg::input_file_validator;

public:
    using base_t::base_t;

    virtual void operator()(std::filesystem::path const & file) const override
    {
        if (file != "-" && file != "/dev/stdin")
            sharg::input_file_validator::operator()(file);
    }
};

class output_file_or_stdout_validator : public sharg::output_file_validator
{
private:
    using base_t = sharg::output_file_validator;

public:
    using base_t::base_t;

    virtual void operator()(std::filesystem::path const & file) const override
    {
        if (file != "-" && file != "/dev/stdout")
            sharg::output_file_validator::operator()(file);
    }
};

// ============================================================================
// Initialisation and program setup
// ============================================================================

// TODO we need to move more of this into bioc++
inline auto create_writer(std::filesystem::path const & filename, char format, size_t const threads)
{
    bool to_stdout = filename == "-" || filename == "/dev/stdout";

    if (to_stdout && format == 'a')
        format = 'v';

    bio::io::var::writer_options writer_opts{.stream_options =
                                               bio::io::transparent_ostream_options{.threads = threads + 1}};

    std::variant<bio::io::bcf, bio::io::vcf> var;

    switch (format)
    {
        case 'a':
            return bio::io::var::writer{filename, writer_opts};
        case 'b':
            var                                    = bio::io::bcf{};
            writer_opts.stream_options.compression = bio::io::compression_format::bgzf;
            break;
        case 'u':
            var                                    = bio::io::bcf{};
            writer_opts.stream_options.compression = bio::io::compression_format::none;
            break;
        case 'z':
            var                                    = bio::io::vcf{};
            writer_opts.stream_options.compression = bio::io::compression_format::bgzf;
            break;
        case 'v':
            var                                    = bio::io::vcf{};
            writer_opts.stream_options.compression = bio::io::compression_format::none;
            break;
        default:
            BIOCPP_UNREACHABLE;
    }

    if (to_stdout)
        return bio::io::var::writer{std::cout, var, writer_opts};
    else
        return bio::io::var::writer{filename, var, writer_opts};
}

template <typename T>
inline void concatenated_sequences_create_scaffold(bio::ranges::concatenated_sequences<T> & concat_seqs,
                                                   size_t const                             outer_size,
                                                   size_t const                             inner_size)
{
    concat_seqs.clear();
    auto && [data_vec, data_delim] = concat_seqs.raw_data();

    data_vec.resize(outer_size * inner_size);

    data_delim.resize(outer_size + 1ul);
    for (size_t i = 0; i < data_delim.size(); ++i)
        data_delim[i] = i * inner_size;

    assert(data_delim.back() == data_vec.size());
}
