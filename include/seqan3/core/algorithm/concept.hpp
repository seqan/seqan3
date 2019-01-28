// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for the configuration classes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{
//!\cond
// Forward declarations
template <typename derived_t, typename value_t>
struct pipeable_config_element;
//!\endcond
}

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Concept config_element_concept
// ----------------------------------------------------------------------------

//!\cond
template <typename derived_t>
class config_element_base;
//!\endcond

/*!\interface seqan3::detail::config_element_concept <>
 * \brief Concept for an algorithm configuration.
 * \ingroup algorithm
 *
 * \extends    std::Semiregular
 * \implements seqan3::pipeable_config_element
 */

/*!\name Requirements for seqan3::detail::config_element_concept
 * \relates seqan3::detail::config_element_concept
 * \brief   You can expect this member on all types that satisfy seqan3::detail::config_element_concept.
 * \{
 */
/*!\var id
 * \memberof seqan3::detail::config_element_concept
 * \brief Algorithm specific static id used for internal validation checks.
 */
/*!\var value
 * \memberof seqan3::detail::config_element_concept
 * \brief Member storing the configuration value.
 */
//!\}
//!\cond
template <typename config_t>
SEQAN3_CONCEPT config_element_concept = std::Semiregular<std::remove_reference_t<config_t>> &&
requires (config_t c)
{
    { std::remove_reference_t<config_t>::id };
    { c.value };
    // Must inherit from the pipeable_config_element class.
    requires std::is_base_of_v<pipeable_config_element<config_t, decltype(c.value)>, config_t>;
};
//!\endcond
} // namespace seqan3::detail
