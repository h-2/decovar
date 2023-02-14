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

// std::vector{std::move(foo), std::move(bar)} actually creates local copies of foo and bar :@
template <typename T, typename... Args>
auto make_vector(T && arg1, Args &&... args)
{
    std::vector<T> vec;
    vec.reserve(sizeof...(Args) + 1);
    vec.emplace_back(std::forward<T>(arg1));
    (vec.emplace_back(std::forward<Args>(args)), ...);
    return vec;
}
