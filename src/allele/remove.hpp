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

namespace _remove
{

/* vectors for fields of multiplicity A, R or G indicating whether that value at the positions shall be removed(1)
* or not (0) */
struct cache_t
{
    std::vector<int /*bool*/> A;
    std::vector<int /*bool*/> R;
    std::vector<int /*bool*/> G;

    std::vector<std::pair<size_t, size_t>> formula_reverse_cache;
};

/* this needs to be run AFTER the filter-vector-R has been computed */
inline void determine_filter_vector_AG(size_t const n_alts, cache_t & filter_vectors) // <- in-out-param
{
    /* filtered alleles A */
    filter_vectors.A.resize(n_alts);
    std::ranges::copy(filter_vectors.R.begin() + 1, filter_vectors.R.end(), filter_vectors.A.begin());

    /* filtered alleles G */
    size_t const gt_size = bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1;
    filter_vectors.G.resize(gt_size);
    for (size_t b = 0; b <= n_alts; ++b)
    {
        for (size_t a = 0; a <= b; ++a)
        {
            assert(bio::io::var::detail::vcf_gt_formula(a, b) < filter_vectors.G.size());
            filter_vectors.G[bio::io::var::detail::vcf_gt_formula(a, b)] = filter_vectors.R[a] || filter_vectors.R[b];
        }
    }

    /* formula_reverse_cache */
    if (filter_vectors.formula_reverse_cache.size() < gt_size)
    {
        filter_vectors.formula_reverse_cache.resize(gt_size);
        for (size_t b = 0; b <= n_alts; ++b)
            for (size_t a = 0; a <= b; ++a)
                filter_vectors.formula_reverse_cache[bio::io::var::detail::vcf_gt_formula(a, b)] = {a, b};
    }
}

inline void determine_filter_vector_R(record_t::info_t const & record_info,
                                      size_t const             record_no,
                                      size_t const             n_alts,
                                      program_options const &  opts,
                                      cache_t &                filter_vectors) // <- out-param
{
    bool has_AF = false;

    /* filtered alleles R */
    filter_vectors.R.clear();
    filter_vectors.R.resize(n_alts + 1);
    filter_vectors.R[0] = 0;
    for (auto && [id, value] : record_info)
    {
        if (id == "AF")
        {
            if (!std::holds_alternative<std::vector<float>>(value))
            {
                throw decovar_error{
                  "[Record no: {}] AF field of multi-allelic record wasn't "
                  "vector<float>.",
                  record_no};
            }

            std::vector<float> const & afs = std::get<std::vector<float>>(value);
            if (afs.size() != n_alts)
            {
                throw decovar_error{
                  "[Record no: {}] AF field of multi-allelic record has wrong "
                  "size: {}, but {} was expected.",
                  record_no,
                  afs.size(),
                  n_alts};
            }

            for (size_t i = 0; i < n_alts; ++i)
                if (afs[i] < opts.rare_af_threshold)
                    filter_vectors.R[i + 1] = 1; // filtered_alleles[0] is REF and never filtered

            has_AF = true;
            break;
        }
    }

    if (!has_AF)
    {
        //TODO look for AC and AN and update those?
        throw decovar_error{"[Record no: {}] no AF field in record.", record_no};
    }
}

template <typename T>
inline void remove_by_indexes(std::vector<T> & vec, std::span<int const> const filter_vector)
{
    auto pred = [&](T const & elem) -> bool
    {
        ptrdiff_t i = (&elem - &*vec.begin()) % filter_vector.size(); // modulo down-maps for concat's inner vector
        assert(i >= 0);
        return filter_vector[i] != 0;
    };
    auto removed_range = std::ranges::remove_if(vec.begin(), vec.end(), pred);
    vec.resize(vec.size() - removed_range.size());
}

inline void update_infos(record_t::info_t & record_info, //??? in-out parameter
                         header_t const &   hdr,
                         size_t const       record_no,
                         cache_t const &    filter_vectors)
{
    std::span<int const> selected_filter_vector{};

    for (auto && [_id, value] : record_info)
    {
        std::string_view         id   = _id;
        size_t const             i    = hdr.string_to_info_pos().at(id);
        header_t::info_t const & info = hdr.infos[i];

        auto visitor = bio::meta::overloaded{
          [&](auto) {
              throw decovar_error{"[Record no: {}] Expected a vector when trimming field {}.", record_no, id};
          },
          [&]<typename T>(std::vector<T> & vec)
          {
              if (vec.size() != selected_filter_vector.size())
              {
                  throw decovar_error{
                    "[Record no: {}] Expected {} elements in field {}, but got {}. A single '.' "
                    "as placeholder is currently not supported.",
                    selected_filter_vector.size(),
                    id,
                    vec.size(),
                    record_no};
              }

              remove_by_indexes(vec, selected_filter_vector);
#ifndef NDEBUG
              size_t ones = std::accumulate(selected_filter_vector.begin(), selected_filter_vector.end(), 0ul);
#endif
              assert(vec.size() == selected_filter_vector.size() - ones);
          }};

        switch (info.number)
        {
            case bio::io::var::header_number::R:
                selected_filter_vector = filter_vectors.R;
                std::visit(visitor, value);
                break;
            case bio::io::var::header_number::A:
                selected_filter_vector = filter_vectors.A;
                std::visit(visitor, value);
                break;
            default:
                break;
        }
    }
}

inline void update_genotypes(record_t::genotypes_t & record_genotypes, //??? in-out parameter
                             header_t const &        hdr,
                             size_t const            record_no,
                             cache_t const &         filter_vectors)
{
    std::span<int const> selected_filter_vector{};

    for (auto && [_id, value] : record_genotypes)
    {
        std::string_view           id     = _id;
        size_t const               i      = hdr.string_to_format_pos().at(id);
        header_t::format_t const & format = hdr.formats[i];

        auto visitor = bio::meta::overloaded{
          [&](auto) {
              throw decovar_error{"[Record no: {}] Expected a vector when trimming field {}.", record_no, id};
          },
          [&]<typename T>(bio::ranges::concatenated_sequences<std::vector<T>> & vec)
          {
              size_t const n_samples        = vec.size();
              size_t const n_alleles_before = selected_filter_vector.size();
              size_t const n_alleles_after =
                n_alleles_before - std::accumulate(selected_filter_vector.begin(), selected_filter_vector.end(), 0ul);

              if (vec.concat_size() != n_samples * n_alleles_before)
              {
                  throw decovar_error{
                    "[Record no: {}] Currently, every sample must be diploid and must contain the "
                    "the correct number of values (e.g. no single '.' placeholder allowed).",
                    record_no};
              }

              std::pair raw_data = vec.raw_data();

              assert(vec.concat_size() == n_samples * n_alleles_before);
              assert(raw_data.second.back() == vec.concat_size());

              remove_by_indexes(raw_data.first, selected_filter_vector);
              assert(vec.concat_size() == n_samples * n_alleles_after);

              assert(vec.size() == n_samples);                 // size remains unchanged
              assert(raw_data.second.size() == n_samples + 1); // size remains unchanged

              for (size_t i = 0; i < n_samples + 1; ++i)
              {
                  assert(raw_data.second[i] == i * n_alleles_before);
                  raw_data.second[i] = i * n_alleles_after;
              }
              assert(raw_data.second.back() == raw_data.first.size());

              if (id == "PL") // PL values are renormalised, so the smallest PL value is 0
              {
                  for (std::span<T> const sample_PL : vec)
                  {
                      T const min = *std::ranges::min_element(sample_PL);
                      if (min > 0)
                          for (T & PL_value : sample_PL)
                              PL_value -= min;
                  }
              }
          }};

        switch (format.number)
        {
            case bio::io::var::header_number::R:
                selected_filter_vector = filter_vectors.R;
                std::visit(visitor, value);
                break;
            case bio::io::var::header_number::A:
                selected_filter_vector = filter_vectors.A;
                std::visit(visitor, value);
                break;
            case bio::io::var::header_number::G:
                selected_filter_vector = filter_vectors.G;
                std::visit(visitor, value);
                break;
            default:
                break;
        }
    }
}

inline void fix_GT(record_t::genotypes_t & record_genotypes, //??? in-out parameter
                   size_t const            record_no,
                   cache_t const &         filter_vectors)
{
    std::span<std::string> all_GT{};

    for (auto && [id, value] : record_genotypes)
        if (id == "GT")
            all_GT = std::get<std::vector<std::string>>(value), ({ break; });

    if (all_GT.empty()) // no GT field present
        return;

    auto visitor =
      bio::meta::overloaded{[&](auto) {
                                throw decovar_error{"[Record no: {}] Expected a vector when reading PL.", record_no};
                            },
                            [&]<typename T>(bio::ranges::concatenated_sequences<std::vector<T>> & vec)
                            {
                                for (size_t i = 0; i < vec.size(); ++i)
                                {
                                    std::span<T const> const sample_PL = vec[i];
                                    std::string &            sample_GT = all_GT[i];

                                    size_t const i_min = std::ranges::min_element(sample_PL) - sample_PL.begin();

                                    auto [a, b] = filter_vectors.formula_reverse_cache[i_min];
                                    sample_GT.clear();
                                    fmt::format_to(std::back_inserter(sample_GT), "{}/{}", a, b); // always unphased
                                }
                            }};

    for (auto && [id, value] : record_genotypes)
        if (id == "PL")
            std::visit(visitor, value), ({ break; });
}

// returns true if all alleles were removed and the entire record should be skipped
[[nodiscard]] inline bool remove_rare_alleles(record_t &              record,
                                              size_t const            record_no,
                                              header_t const &        hdr,
                                              program_options const & opts,
                                              cache_t &               filter_vectors)
{
    size_t n_alts = record.alt.size();

    determine_filter_vector_R(record.info, record_no, n_alts, opts, filter_vectors);
    determine_filter_vector_AG(n_alts, filter_vectors);

    log(opts, "filter_vector.A: {}\n", filter_vectors.A);
    log(opts, "filter_vector.R: {}\n", filter_vectors.R);
    log(opts, "filter_vector.G: {}\n", filter_vectors.G);

    if (std::ranges::all_of(filter_vectors.A, std::identity{})) // all alleles are removed ??? skip this record
    {
        log(opts, "record no {} would have no remaining alleles and is skipped completely.\n", record_no);
        return true;
    }

    if (std::ranges::any_of(filter_vectors.A, std::identity{})) // only modify record if necessary
    {
        /* update alts */
        remove_by_indexes(record.alt, filter_vectors.A);

        /* update info */
        update_infos(record.info, hdr, record_no, filter_vectors);

        /* update genotypes */
        update_genotypes(record.genotypes, hdr, record_no, filter_vectors);

        /* fix GT values after alleles have been removed */
        fix_GT(record.genotypes, record_no, filter_vectors);
    }

    return false;
}

} // namespace _remove
