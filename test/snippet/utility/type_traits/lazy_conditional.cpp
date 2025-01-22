// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

//![complete]
#include <forward_list> // std::forward_list
#include <ranges>       // std::ranges::input_range
#include <vector>       // std::vector

#include <seqan3/core/range/type_traits.hpp>               // seqan3::size_type_t
#include <seqan3/utility/type_traits/lazy_conditional.hpp> // seqan3::lazy_conditional_t

template <std::ranges::input_range rng_t>
void foobar(rng_t && range)
{

#if 0
    // The following would fail to compile if rng_t is not sized,
    // because std::ranges::range_size_t<rngt_t> needs to be valid
    // (independent of whether the condition evaluates to true)
    using size_type = std::conditional_t<std::ranges::sized_range<rng_t>,
                                         std::ranges::range_size_t<rng_t>,
                                         void>;
#endif

    // This delays instantiation of std::ranges::range_size_t<rngt_t> until after the
    // conditional-decision is made:
    using size_type = seqan3::detail::lazy_conditional_t<std::ranges::sized_range<rng_t>,
                                                         seqan3::detail::lazy<std::ranges::range_size_t, rng_t>,
                                                         void>;

    // DO SOMETHING with size_type
}

int main()
{
    foobar(std::vector<int>{});       // sized
    foobar(std::forward_list<int>{}); // not sized
}
//![complete]

#pragma GCC diagnostic pop
