// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for the configuration classes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{
//!\cond
// Forward declarations
template <typename derived_t, typename value_t = void>
struct pipeable_config_element;
//!\endcond
}

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Concept config_element_specialisation
// ----------------------------------------------------------------------------

//!\cond
template <typename derived_t>
class config_element_base;
//!\endcond

/*!\interface seqan3::detail::config_element_specialisation <>
 * \brief Concept for an algorithm configuration.
 * \ingroup algorithm
 *
 * \extends    std::semiregular
 * \implements seqan3::pipeable_config_element
 */

/*!\name Requirements for seqan3::detail::config_element_specialisation
 * \relates seqan3::detail::config_element_specialisation
 * \brief   You can expect this member on all types that satisfy seqan3::detail::config_element_specialisation.
 * \{
 */
/*!\var id
 * \brief Algorithm specific static id used for internal validation checks.
 */
//!\}
//!\cond
template <typename config_t>
SEQAN3_CONCEPT config_element_specialisation = std::semiregular<std::remove_reference_t<config_t>> &&
requires (config_t c)
{
    { std::remove_reference_t<config_t>::id };
};
//!\endcond

} // namespace seqan3::detail
