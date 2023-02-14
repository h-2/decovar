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

#include <ranges>
#include <variant>

#include <bio/ranges/container/concatenated_sequences.hpp>

#include "../misc.hpp"
#include "allele.hpp"
#include "remove.hpp"

namespace _split
{

enum which_t
{
    leq,
    gt
};

inline void determine_filter_vector_R(record_t const &        record,
                                      program_options const & opts,
                                      which_t const           which,     // if true, remove lengths LEQ option
                                      _remove::cache_t &      filter_vectors) // <- out-param
{
    size_t const n_alts = record.alt.size();

    /* filtered alleles R */
    filter_vectors.R.clear();
    filter_vectors.R.resize(n_alts + 1);
    filter_vectors.R[0] = 0; // REF is never filtered
    for (size_t i = 0; i < n_alts; ++i)
        filter_vectors.R[i + 1] = (record.alt[i].size() <= opts.split_by_length) == (which == leq);
}

bool needs_splitting(record_t const & record, program_options const & opts)
{
    size_t n_alts = record.alt.size();
    if (n_alts <= 1)
        return false;

    bool has_shorter = false;
    bool has_longer  = false;

    for (auto const & alt_allele : record.alt)
    {
        if (alt_allele.size() <= opts.split_by_length)
            has_shorter = true;
        if (alt_allele.size() > opts.split_by_length)
            has_longer = true;
    }

    return has_shorter && has_longer;
}

// returns true if all alleles were removed and the entire record should be skipped
inline void remove_alleles(record_t &              record,
                           size_t const            record_no,
                           which_t const           which,
                           header_t const &        hdr,
                           program_options const & opts,
                           _remove::cache_t &      filter_vectors)
{
    size_t n_alts = record.alt.size();

    determine_filter_vector_R(record, opts, which, filter_vectors);
    _remove::determine_filter_vector_AG(n_alts, filter_vectors);

    log(opts, "filter_vector.A: {}\n", filter_vectors.A);
    log(opts, "filter_vector.R: {}\n", filter_vectors.R);
    log(opts, "filter_vector.G: {}\n", filter_vectors.G);

    /* update alts */
    _remove::remove_by_indexes(record.alt, filter_vectors.A);

    /* update info */
    _remove::update_infos(record.info, hdr, record_no, filter_vectors);

    /* update genotypes */
    _remove::update_genotypes(record.genotypes, hdr, record_no, filter_vectors);

    /* fix GT values after alleles have been removed */
    _remove::fix_GT(record.genotypes, record_no, filter_vectors);
}

} // namespace _split
