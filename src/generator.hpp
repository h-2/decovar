#pragma once

#include <version>

#if defined(__cpp_lib_generator) && __has_include(<generator>)
#    include <generator>
#else
#    include "../submodules/generator/include/__generator.hpp"
#endif

#include <bio/ranges/views/detail.hpp>

inline constexpr auto cojoin =
  []<std::ranges::input_range rng_t>(rng_t &&
                                     urange) requires(std::ranges::input_range<std::ranges::range_reference_t<rng_t>>)
{
    return [](auto v) -> std::generator<std::ranges::range_reference_t<std::ranges::range_reference_t<rng_t>>>
    {
        for (auto && subrange : v)
            co_yield std::ranges::elements_of(std::forward<decltype(subrange)>(subrange));
    }(std::views::all(std::forward<rng_t>(urange)));
};

inline constexpr auto views_cojoin = bio::ranges::detail::adaptor_from_functor{cojoin};
