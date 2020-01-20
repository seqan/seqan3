// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides my_tuple for testing tuple utility.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <utility>

#include <seqan3/core/type_traits/basic.hpp>

namespace seqan3
{
struct my_tuple
{
    int   el0;
    float el1;

    constexpr bool operator==(my_tuple const & rhs) const
    {
        return std::tie(el0, el1) == std::tie(rhs.el0, rhs.el1);
    }

    constexpr bool operator!=(my_tuple const & rhs) const
    {
        return !(*this == rhs);
    }

    constexpr bool operator<(my_tuple const & rhs) const
    {
        return std::tie(el0, el1) < std::tie(rhs.el0, rhs.el1);
    }

    constexpr bool operator<=(my_tuple const & rhs) const
    {
        return std::tie(el0, el1) <= std::tie(rhs.el0, rhs.el1);
    }

    constexpr bool operator>(my_tuple const & rhs) const
    {
        return std::tie(el0, el1) > std::tie(rhs.el0, rhs.el1);
    }

    constexpr bool operator>=(my_tuple const & rhs) const
    {
        return std::tie(el0, el1) >= std::tie(rhs.el0, rhs.el1);
    }
};

template <size_t elem>
constexpr auto & get(seqan3::my_tuple & t)
{
    static_assert(elem < 2);

    if constexpr (elem == 0)
        return t.el0;
    else
        return t.el1;
}

template <size_t elem>
constexpr auto const & get(seqan3::my_tuple const & t)
{
    static_assert(elem < 2);

    if constexpr (elem == 0)
        return t.el0;
    else
        return t.el1;
}

template <size_t elem>
constexpr auto && get(seqan3::my_tuple && t)
{
    static_assert(elem < 2);

    if constexpr (elem == 0)
        return std::move(t.el0);
    else
        return std::move(t.el1);
}

template <size_t elem>
constexpr auto const && get(seqan3::my_tuple const && t)
{
    static_assert(elem < 2);

    if constexpr (elem == 0)
        return std::move(t.el0);
    else
        return std::move(t.el1);
}

} // namespace seqan3

namespace std
{

template <>
struct tuple_size<seqan3::my_tuple>
{
    static constexpr size_t value = 2;
};

template <size_t elem_no>
struct tuple_element<elem_no, seqan3::my_tuple>
{
    using type = seqan3::remove_cvref_t<decltype(get<elem_no>(std::declval<seqan3::my_tuple>()))>;
};

} // namespace std
