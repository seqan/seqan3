// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <type_traits>

#include <seqan3/utility/type_pack/detail/type_pack_algorithm.hpp>

// With c++20 you could also write it like this
// auto fn = []<typename value_t>(value_t && value)
// {
// ...
// };
auto fn = [](auto value)
{
    // id is the original type not wrapped in std::type_identity.
    using value_t = decltype(value);

    if constexpr (std::is_same_v<value_t, bool>)
        return value == false;
    else if constexpr (std::is_same_v<value_t, int>)
        return value == 3;
    else if constexpr (std::is_same_v<value_t, double>)
        return value - 1.2 < 0.00001;
    else
        return false;
};

static_assert(seqan3::detail::all_of(fn, 3, 1.2, false));                    // evaluates to true
static_assert(!seqan3::detail::all_of(fn, 3, 1.2, false, "something else")); // evaluates to false
