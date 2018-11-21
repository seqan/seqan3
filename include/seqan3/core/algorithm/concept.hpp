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
 * \brief Provides concepts for the configuration classes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Concept ConfigElement
// ----------------------------------------------------------------------------

//!\cond
template <typename derived_t>
class config_element_base;
//!\endcond

/*!\interface seqan3::detail::ConfigElement <>
 * \brief Concept for an algorithm configuration.
 * \ingroup algorithm
 *
 * \extends std::Semiregular
 */

/*!\name Requirements for seqan3::detail::ConfigElement
 * \relates seqan3::detail::ConfigElement
 * \brief   You can expect this member on all types that satisfy seqan3::detail::ConfigElement.
 * \{
 */
/*!\var value
 * \memberof seqan3::detail::ConfigElement
 * \brief Member storing the configuration value.
 */
//!\}
//!\cond
template <typename config_t>
concept ConfigElement = std::Semiregular<std::remove_reference_t<config_t>> &&
requires (config_t c)
{
    { c.value };
};
//!\endcond

// ----------------------------------------------------------------------------
// Concept DeferredConfigElement
// ----------------------------------------------------------------------------

//!\cond
// Forward declaration for the deferred config concept
template <detail::ConfigElement ... configs_t>
class configuration;

template <typename derived_t>
class deferred_config_element_base;
//!\endcond

/*!\interface seqan3::detail::DeferredConfigElement <>
 * \brief Concept for a deferred algorithm configuration.
 * \ingroup algorithm
 *
 * \extends seqan3::detail::ConfigElement
 * \implements seqan3::detail::deferred_config_element_base
 *
 * \details
 *
 * Classes must extend the abstract crtp-base class seqan3::detail::deferred_config_element_base in order to satisfy
 * this concept. This concept is merely used for internal purposes and shall not be exposed in public user interfaces.
 */
//!\cond
template <typename config_t>
concept DeferredConfigElement = ConfigElement<config_t> &&
    std::is_base_of_v<deferred_config_element_base<std::remove_reference_t<config_t>>,
                                                   std::remove_reference_t<config_t>>;
//!\endcond

// ----------------------------------------------------------------------------
// Metafunction is_algorithm_configuration
// ----------------------------------------------------------------------------

/*!\brief Value metafunction that returns whether a type is an algorithm configuration.
 * \ingroup algorithm
 *
 * \returns std::true_type if the given type is a seqan3::detail::configuration, else std::false_type.
 */
template <typename object_t>
struct is_algorithm_configuration : std::false_type
{};

//!\cond
template <typename ...config_elements_t>
struct is_algorithm_configuration<configuration<config_elements_t...>> : std::true_type
{};
//!\endcond

/*!\brief Helper variable template for seqan3::detail::is_algorithm_configuration.
 * \ingroup algorithm
 */
template <typename object_t>
inline constexpr bool is_algorithm_configuration_v = is_algorithm_configuration<object_t>::value;

} // namespace seqan3::detail
