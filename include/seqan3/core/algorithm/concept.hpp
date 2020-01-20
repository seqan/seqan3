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
template <typename derived_t, typename value_t>
struct pipeable_config_element;
//!\endcond
}

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Concept config_element
// ----------------------------------------------------------------------------

//!\cond
template <typename derived_t>
class config_element_base;
//!\endcond

/*!\interface seqan3::detail::config_element <>
 * \brief Concept for an algorithm configuration.
 * \ingroup algorithm
 *
 * \extends    std::semiregular
 * \implements seqan3::pipeable_config_element
 */

/*!\name Requirements for seqan3::detail::config_element
 * \relates seqan3::detail::config_element
 * \brief   You can expect this member on all types that satisfy seqan3::detail::config_element.
 * \{
 */
/*!\var id
 * \brief Algorithm specific static id used for internal validation checks.
 */
/*!\var value
 * \brief Member storing the configuration value.
 */
//!\}
//!\cond
template <typename config_t>
SEQAN3_CONCEPT config_element = std::semiregular<std::remove_reference_t<config_t>> &&
requires (config_t c)
{
    { std::remove_reference_t<config_t>::id };
    { c.value };
    // Must inherit from the pipeable_config_element class.
    requires std::is_base_of_v<pipeable_config_element<config_t, decltype(c.value)>, config_t>;
};
//!\endcond
} // namespace seqan3::detail
