// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides my_tuple for testing tuple utility.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <utility>

#include <seqan3/core/metafunction/basic.hpp>

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
