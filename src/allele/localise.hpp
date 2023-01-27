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

#include <cstdint>
#include <ranges>

#include <bio/ranges/container/concatenated_sequences.hpp>
#include <bio/ranges/views/zip.hpp>

#include <bio/io/var/record.hpp>

#include "allele.hpp"
#include "../misc.hpp"

struct localise_cache_t
{
    // vector of size L that contains the indexes of the retained alleles for each sample
    bio::ranges::concatenated_sequences<std::vector<int32_t /*TODO?*/>> laa;

    bio::ranges::concatenated_sequences<std::vector<int8_t>>  vec8;
    bio::ranges::concatenated_sequences<std::vector<int16_t>> vec16;
    bio::ranges::concatenated_sequences<std::vector<int32_t>> vec32;

    template <std::signed_integral int_t>
    auto & get_buf()
    {
        if constexpr (BIOCPP_IS_SAME(int_t, int8_t))
            return vec8;
        else if constexpr (BIOCPP_IS_SAME(int_t, int16_t))
            return vec16;
        else
            return vec32;
    }

    std::vector<std::pair<double, size_t>> probs_buf;

    std::vector<std::pair<int8_t, size_t>>  pair_buf8;
    std::vector<std::pair<int16_t, size_t>> pair_buf16;
    std::vector<std::pair<int32_t, size_t>> pair_buf32;

    template <std::signed_integral int_t>
    auto & get_pair_buf()
    {
        if constexpr (BIOCPP_IS_SAME(int_t, int8_t))
            return pair_buf8;
        else if constexpr (BIOCPP_IS_SAME(int_t, int16_t))
            return pair_buf16;
        else
            return pair_buf32;
    }
};

inline double PL_to_prob(int32_t const PL_val)
{
    return std::pow(10.0, static_cast<double>(PL_val) / -10.0);
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

template <typename int_t>
inline void determine_laa(localise_cache_t &                                              cache,
                          bio::ranges::concatenated_sequences<std::vector<int_t>> const & PLs,
                          record_t const &                                                record,
                          size_t const                                                    record_no,
                          header_t const &                                                hdr,
                          program_options const &                                         opts)
{
    auto &       laa       = cache.laa;
    size_t       n_alts    = record.alt.size();
    size_t       n_samples = hdr.column_labels.size() - 9;
    size_t const L         = opts.local_alleles;

    if (PLs.concat_size() != n_samples * (bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1))
    {
        throw decovar_error{
          "[Record no: {}] Currently, every sample must be diploid and must contain the "
          "full number of PL values (e.g. no single '.' placeholder allowed).",
          record_no};
    }

    laa.clear();
    laa.concat_reserve(n_samples * L);

    for (std::span<int_t const> const sample_PLs : PLs)
    {
        cache.probs_buf.clear();
        cache.probs_buf.resize(n_alts + 1);

        for (size_t i = 0; i < cache.probs_buf.size(); ++i)
            cache.probs_buf[i].second = i;

        for (size_t b = 0; b <= n_alts; ++b)
        {
            for (size_t a = 0; a <= b; ++a)
            {
                assert(bio::io::var::detail::vcf_gt_formula(a, b) < sample_PLs.size());
                double prob = PL_to_prob(sample_PLs[bio::io::var::detail::vcf_gt_formula(a, b)]);

                cache.probs_buf[a].first += prob;
                cache.probs_buf[b].first += prob;
            }
        }

        // we sort all (except REF allele) by probability
        std::ranges::sort(cache.probs_buf.begin() + 1,
                          cache.probs_buf.end(),
                          std::ranges::greater{},
                          [](auto && pair) { return pair.first; });

        // now we sort the first L + 1 by their index again
        std::ranges::sort(cache.probs_buf.begin(),
                          cache.probs_buf.begin() + L + 1,
                          std::ranges::less{},
                          [](auto && pair) { return pair.second; }); // it should be possible to do this nicer

        // 0 (the REF position) is not copied to much dismay; only the next L
        laa.push_back(cache.probs_buf | std::views::elements<1> | bio::views::slice(1, L + 1));
        //NOTE, this is likely not the most efficient solution
    }

    assert(laa.concat_size() == n_samples * L);

    log(opts, "Index map: {}\n", laa);
}

// returns true if all alleles were removed and the entire record should be skipped
inline void localise_alleles(record_t &              record,
                             size_t const            record_no,
                             header_t const &        hdr,
                             program_options const & opts,
                             localise_cache_t &      cache)
{
    size_t const n_alts    = record.alt.size();
    size_t const n_samples = hdr.column_labels.size() - 9;
    size_t const L         = opts.local_alleles;
    assert(n_alts > L);

    // TODO the following can be replaced once we have bio::ranges::dictionary
    std::unordered_map<std::string_view, size_t> field_ids;
    for (size_t i = 0; i < record.genotypes.size(); ++i)
        field_ids[record.genotypes[i].id] = i;

    for (std::string_view id : {"LAA", "LAD", "LGT", "LPL"})
        if (field_ids.contains(id))
            throw decovar_error{"[Record no: {}] Cannot add {} field, because {} field already present.",
                                record_no,
                                id,
                                id};

    if (auto it = field_ids.find("PL"); it != field_ids.end())
    {
        auto fn = bio::meta::overloaded{
          [&]<std::integral int_t>(bio::ranges::concatenated_sequences<std::vector<int_t>> const & PLs)
          { determine_laa(cache, PLs, record, record_no, hdr, opts); },
          [record_no](auto const &) {
              throw decovar_error{"[Record no: {}] PL-field was in wrong state.", record_no};
          }};

        std::visit(fn, record.genotypes[it->second].value);
    }
    else
    {
        throw decovar_error{"[Record no: {}] Cannot compute localised alleles if PL-field is not present.", record_no};
    }

    /* LAD */
    if (auto it = field_ids.find("AD"); it != field_ids.end())
    {
        auto visitor = bio::meta::overloaded{
          [&]<std::integral int_t>(bio::ranges::concatenated_sequences<std::vector<int_t>> & field_AD)
          {
              auto & buffer = cache.get_buf<int_t>();
              buffer.clear();
              buffer.reserve(n_samples);
              buffer.concat_reserve(n_samples * (L + 1));

              assert(field_AD.size() == cache.laa.size());
              for (auto const && [sample_AD, sample_laa] : bio::views::zip(field_AD, cache.laa))
              {
                  buffer.push_back();                   // append empty
                  buffer.push_back_inner(sample_AD[0]); // reference is always appended
                  for (int32_t i : sample_laa)
                      buffer.push_back_inner(sample_AD[i]);
              }
              // NOTE: instead of reserving and appending, resize and assign might be faster

              record.genotypes.emplace_back("LAD", std::move(buffer)); // create LAD field

              if (opts.remove_global_alleles)
                  buffer = std::move(field_AD); // salvage dynamic memory from field_AD since it will be removed later
          },
          [record_no](auto &) {
              throw decovar_error{"[Record no: {}] LAD field was not a range of integers.", record_no};
          }};

        std::visit(visitor, record.genotypes[it->second].value);
    }

    /* LGT */
    // TODO what is the point of this?

    /* LPL */
    if (auto it = field_ids.find("PL"); it != field_ids.end())
    {
        auto visitor = bio::meta::overloaded{
          [&]<std::integral int_t>(bio::ranges::concatenated_sequences<std::vector<int_t>> & field_PL)
          {
              // correct size of PL was already checked above

              auto & buffer = cache.get_buf<int_t>();
              concatenated_sequences_create_scaffold(buffer, n_samples, bio::io::var::detail::vcf_gt_formula(L, L) + 1);

              for (size_t i = 0; i < n_samples; ++i)
              {
                  std::span<int32_t const> sample_LAA = cache.laa[i];
                  std::span<int_t const>   sample_PL  = field_PL[i];
                  std::span<int_t>         sample_LPL = buffer[i];

                  assert(sample_PL.size() == bio::io::var::detail::vcf_gt_formula(n_alts, n_alts) + 1);
                  assert(sample_LPL.size() == bio::io::var::detail::vcf_gt_formula(L, L) + 1);

                  /* NOTE normally we would loop over [0, L]
                   * however, the sample_LAA does not contain 0 at position 0 (for the ref), it only has the
                   * alternative alleles. So that's why we need to subtract 1 in the mapping and why we need
                   * to create special cases for a=0 and/or b=0
                   */
                  sample_LPL[0] = sample_PL[0]; // formula(0,0) is 0 and REF is always preserved
                  for (size_t b = 1; b <= L; ++b)
                  {
                      sample_LPL[bio::io::var::detail::vcf_gt_formula(0, b)] =
                        sample_PL[bio::io::var::detail::vcf_gt_formula(0, sample_LAA[b - 1])];

                      for (size_t a = 1; a <= b; ++a)
                      {
                          assert(bio::io::var::detail::vcf_gt_formula(a, b) < sample_LPL.size());
                          assert(bio::io::var::detail::vcf_gt_formula(sample_LAA[a - 1], sample_LAA[b - 1]) <
                                 sample_PL.size());
                          sample_LPL[bio::io::var::detail::vcf_gt_formula(a, b)] =
                            sample_PL[bio::io::var::detail::vcf_gt_formula(sample_LAA[a - 1], sample_LAA[b - 1])];
                      }
                  }
              }

              record.genotypes.emplace_back("LPL", std::move(buffer)); // create LPL field

              if (opts.remove_global_alleles)
                  buffer = std::move(field_PL); // salvage dynamic memory from field_AD since it will be removed later
          },
          [record_no](auto &) {
              throw decovar_error{"[Record no: {}] LPL field was not a range of integers.", record_no};
          }};

        std::visit(visitor, record.genotypes[it->second].value);
    }

    /* LAA */
    record.genotypes.emplace_back("LAA", std::move(cache.laa)); // this comes last, because cache.laa is used before

    /* remove AD, GT, PL */
    if (opts.remove_global_alleles)
    {
        std::erase_if(record.genotypes,
                      [](decltype(record.genotypes[0]) genotype)
                      { return genotype.id == "AD" /*|| genotype.id == "GT"*/ || genotype.id == "PL"; });
    }
}

void salvage_localise_cache(record_t & record, localise_cache_t & cache)
{
    // LPL and LAD have already been swapped with PL and AD, so they don't need to be salvaged

    for (auto && [id, value] : record.genotypes)
    {
        if (id == "LAA")
        {
            auto & laa = std::get<bio::ranges::concatenated_sequences<std::vector<int32_t>>>(value);
            cache.laa  = std::move(laa);
            break;
        }
    }
}
