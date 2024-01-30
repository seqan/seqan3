// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides alphabet helper concepts.
 */

#pragma once

#include <concepts>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
/*!\interface   seqan3::detail::weakly_equality_comparable_with <>
 * \ingroup alphabet
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \brief       Requires the two operands to be comparable with `==` and `!=` in both directions.
 */
//!\cond
template <class T, class U>
concept weakly_equality_comparable_with =
    requires (std::remove_reference_t<T> const & t, std::remove_reference_t<U> const & u) {
        requires std::convertible_to<decltype(t == u), bool>;
        requires std::convertible_to<decltype(t != u), bool>;
        requires std::convertible_to<decltype(u == t), bool>;
        requires std::convertible_to<decltype(u != t), bool>;
    };
//!\endcond

/*!\interface   seqan3::detail::weakly_ordered_with <>
 * \ingroup alphabet
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \brief       Requires the two operands to be comparable with `<`, `<=`, `>` and `>=` in both directions.
 */
//!\cond
template <typename t1, typename t2>
concept weakly_ordered_with =
    requires (std::remove_reference_t<t1> const & v1, std::remove_reference_t<t2> const & v2) {
        requires std::convertible_to<decltype(v1 < v2), bool>;
        requires std::convertible_to<decltype(v1 <= v2), bool>;
        requires std::convertible_to<decltype(v1 > v2), bool>;
        requires std::convertible_to<decltype(v1 >= v2), bool>;

        requires std::convertible_to<decltype(v2 < v1), bool>;
        requires std::convertible_to<decltype(v2 <= v1), bool>;
        requires std::convertible_to<decltype(v2 > v1), bool>;
        requires std::convertible_to<decltype(v2 >= v1), bool>;
    };
//!\endcond

} // namespace seqan3::detail
