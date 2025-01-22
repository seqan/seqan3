// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the concepts seqan3::transformation_trait and seqan3::unary_type_trait.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\interface seqan3::transformation_trait
 * \ingroup utility_type_traits
 * \brief Concept for a transformation trait.
 *
 * An object is a transformation trait if it exposes a member type called `type`.
 */
//!\cond
template <typename t>
concept transformation_trait = requires { typename t::type; };
//!\endcond

/*!\interface seqan3::unary_type_trait
 * \ingroup utility_type_traits
 * \brief Concept for a unary traits type.
 *
 * An object is a unary traits type if it is derived from std::integral_constant.
 */
//!\cond
template <typename t>
concept unary_type_trait = std::is_base_of_v<std::integral_constant<typename t::value_type, t::value>, t>;
//!\endcond

} // namespace seqan3
