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
 * \brief Define various helper functionality.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm> // std::lexicographical_compare

#include <seqan3/range/concept.hpp>
#include <seqan3/core/concept/core.hpp>

#include <range/v3/algorithm/equal.hpp> // ranges::equal

namespace seqan3
{

/*!\name Generic comparison operators for two ranges of the same type (ignoring const).
 *\{
 * \tparam lhs_range_type The first input range to compare (seqan3::input_range_concept).
 * \tparam rhs_range_type The second input range to compare (seqan3::input_range_concept).
 *
 * The two input ranges must satisfy the seqan3::input_range_concept and are
 * required to be of the same type, excluding difference in constness.
 *
 * Note: The operator==/!= delegates to the range_v3::equal function and
 *       requires the underlying value type to be weakly comparable
 *       (seqan3::weakly_equality_comparable_concept).
 * Note: The operator <,<=,>,>= delegates to std::lexicographical_compare and
 *       requires the underlying value type to be strictly less comparable (a < b).
 *
 * ### Exception
 *
 * Strong exception is guaranteed since the input ranges are passed by const
 * and no further objects are created.
 *
 * ### Complexity
 *
 * The complexity is linear in the number of elements in the input range but
 * must be multiplied with the complexity of the comparator of the underlying
 * value_type. (e.g. O(n) for std::vector, O(n^2) for std::vector<std::vector>)
 *
 * ### Thread Safety
 *
 * Since all input parameters are passed by const reference the function is
 * considered thread safe.
 *
 * \sa [ranges::equal] (http://ericniebler.github.io/range-v3/group__group-views.html#gafde1fc2968ed1c0a30146a7231bb64fc)
 * \sa [std::lexicographical_compare] (http://en.cppreference.com/w/cpp/algorithm/lexicographical_compare)
 */
//!\brief Returns `true` if \p lhs == \p rhs, `false` otherwise.
template<input_range_concept lhs_range_type, input_range_concept rhs_range_type>
//!\cond
    requires std::is_same_v<std::remove_const<lhs_range_type>, std::remove_const<rhs_range_type>> &&
    weakly_equality_comparable_concept<typename lhs_range_type::value_type, typename rhs_range_type::value_type>
//!\endcond
constexpr bool operator==(lhs_range_type const & lhs, rhs_range_type const & rhs)
{
    return ranges::equal(lhs, rhs);
}

//!\brief Returns `true` if \p lhs != \p rhs, `false` otherwise.
template<input_range_concept lhs_range_type, input_range_concept rhs_range_type>
//!\cond
    requires std::is_same_v<std::remove_const<lhs_range_type>, std::remove_const<rhs_range_type>> &&
    weakly_equality_comparable_concept<typename lhs_range_type::value_type, typename rhs_range_type::value_type>
//!\endcond
constexpr bool operator!=(lhs_range_type const & lhs, rhs_range_type const & rhs)
{
    return !(lhs == rhs);
}

//!\brief Returns `true` if \p lhs < \p rhs, `false` otherwise.
template<input_range_concept lhs_range_type, input_range_concept rhs_range_type>
//!\cond
    requires std::is_same_v<std::remove_const<lhs_range_type>, std::remove_const<rhs_range_type>> &&
             requires (typename lhs_range_type::value_type lv, typename rhs_range_type::value_type rv)
                 { { lv < rv} -> bool; }
//!\endcond
constexpr bool operator<(lhs_range_type const & lhs, rhs_range_type const & rhs)
{
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

//!\brief Returns `true` if \p lhs > \p rhs, `false` otherwise.
template<input_range_concept lhs_range_type, input_range_concept rhs_range_type>
//!\cond
    requires std::is_same_v<std::remove_const<lhs_range_type>, std::remove_const<rhs_range_type>> &&
             requires (typename lhs_range_type::value_type lv, typename rhs_range_type::value_type rv)
                 { { lv < rv} -> bool; }
//!\endcond
constexpr bool operator>(lhs_range_type const & lhs, rhs_range_type const & rhs)
{
    return std::lexicographical_compare(rhs.begin(), rhs.end(), lhs.begin(), lhs.end());
}

//!\brief Returns `true` if \p lhs <= \p rhs, `false` otherwise.
template<input_range_concept lhs_range_type, input_range_concept rhs_range_type>
//!\cond
    requires std::is_same_v<std::remove_const<lhs_range_type>, std::remove_const<rhs_range_type>> &&
             requires (typename lhs_range_type::value_type lv, typename rhs_range_type::value_type rv)
                 { { lv < rv} -> bool; }
//!\endcond
constexpr bool operator<=(lhs_range_type const & lhs, rhs_range_type const & rhs)
{
    return !(lhs > rhs);
}

//!\brief Returns `true` if \p lhs >= \p rhs, `false` otherwise.
template<input_range_concept lhs_range_type, input_range_concept rhs_range_type>
//!\cond
    requires std::is_same_v<std::remove_const<lhs_range_type>, std::remove_const<rhs_range_type>> &&
             requires (typename lhs_range_type::value_type lv, typename rhs_range_type::value_type rv)
                 { { lv < rv} -> bool; }
//!\endcond
constexpr bool operator>=(lhs_range_type const & lhs, rhs_range_type const & rhs)
{
    return !(lhs < rhs);
}
//!\}

} // namespace seqan3
