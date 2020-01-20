// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the concepts seqan3::transformation_trait and seqan3::unary_type_trait.
 * \author Svenja Mehringer <avenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\interface seqan3::transformation_trait
 * \ingroup type_traits
 * \brief Concept for a transformation trait.
 *
 * An object is a transformation trait if it exposes a member type called `type`.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT transformation_trait = requires { typename t::type; };
//!\endcond

/*!\interface seqan3::unary_type_trait
 * \ingroup type_traits
 * \brief Concept for a unary traits type.
 *
 * An object is a unary traits type if it is derived from std::integral_constant.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT unary_type_trait = std::is_base_of_v<std::integral_constant<typename t::value_type, t::value>, t>;
//!\endcond

} // namespace seqan3
