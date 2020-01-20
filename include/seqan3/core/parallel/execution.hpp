// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides execution policies.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

/*!\defgroup execution_policy execution
 * \ingroup parallel
 * \brief Adopts the [STL execution policies](https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t).
 */

namespace seqan3
{
/*!\addtogroup execution_policy
 * \{
 */
/*!\brief Sequenced execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t
 */
class sequenced_policy
{};

/*!\brief Parallel execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t
 */
class parallel_policy
{};

/*!\brief Parallel and unsequenced execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t
 */
class parallel_unsequenced_policy
{};

/*!\brief Unsequenced execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t
 */
class unsequenced_policy
{};

/*!\brief Global execution policy object for sequenced execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag
 */
inline constexpr sequenced_policy            seq{};

/*!\brief Global execution policy object for parallel execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag
 */
inline constexpr parallel_policy             par{};

/*!\brief Global execution policy object for parallel and unsequenced execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag
 */
inline constexpr parallel_unsequenced_policy par_unseq{};

/*!\brief Global execution policy object for unsequenced execution policy.
 * \sa https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag
 */
inline constexpr unsequenced_policy          unseq{};

/*!\brief Test whether a class represents an execution policy.
 * \see seqan3::is_execution_policy_v
 */
template <class t>
struct is_execution_policy : std::false_type
{};

//!\cond
template <>
struct is_execution_policy<sequenced_policy> : std::true_type
{};

template <>
struct is_execution_policy<parallel_policy> : std::true_type
{};

template <>
struct is_execution_policy<parallel_unsequenced_policy> : std::true_type
{};

template <>
struct is_execution_policy<unsequenced_policy> : std::true_type
{};
//!\endcond

/*!\brief Helper variable template for std::is_execution_policy.
 * \relates seqan3::is_execution_policy
 */
template <class t>
inline constexpr bool is_execution_policy_v = is_execution_policy<t>::value;
//!\}
} // namespace seqan3
