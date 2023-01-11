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

#include <sharg/all.hpp>
#include <thread>

#include <bio/io/var/header.hpp>
#include <bio/io/var/record.hpp>

using record_t = bio::io::var::record_default;
using header_t = bio::io::var::header;

struct decovar_error : std::runtime_error
{
    template <typename... T>
    decovar_error(fmt::format_string<T...> format, T &&... args) :
      std::runtime_error{fmt::format(format, std::forward<T>(args)...)}
    {}
};

struct program_options
{
    std::string input_file;
    std::string output_file = "-";

    float  rare_af_threshold     = 0ul;
    size_t local_alleles         = 0ul;
    bool   remove_global_alleles = true;

    size_t threads = std::max<size_t>(2, std::min<size_t>(8, std::thread::hardware_concurrency()));

    bool verbose = false;
};

template <typename... T>
inline void log(program_options const & opts, fmt::format_string<T...> format, T &&... args)
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
