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
concept config_element_concept = std::Semiregular<std::remove_reference_t<config_t>> &&
requires (config_t c)
{
    { std::remove_reference_t<config_t>::id };
    { c.value };
    // Must inherit from the pipeable_config_element class.
    requires std::is_base_of_v<pipeable_config_element<config_t, decltype(c.value)>, config_t>;
};
//!\endcond
} // namespace seqan3::detail


namespace seqan3
{
//!\cond
// Forward declaration for concept definition.
template <detail::config_element_concept ... configs_t>
class configuration;
//!\endcond
} // namespace seqan3

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// Concept deferred_config_element_concept
// ----------------------------------------------------------------------------

//!\cond
template <typename derived_t>
class deferred_config_element_base;

template <detail::config_element_concept ... configs_t>
class configuration;
//!\endcond

/*!\interface seqan3::detail::deferred_config_element_concept <>
 * \brief Concept for a deferred algorithm configuration.
 * \ingroup algorithm
 *
 * \extends seqan3::detail::config_element_concept
 * \implements seqan3::detail::deferred_config_element_base
 *
 * \details
 *
 * Classes must extend the abstract crtp-base class seqan3::detail::deferred_config_element_base in order to satisfy
 * this concept. This concept is merely used for internal purposes and shall not be exposed in public user interfaces.
 */
//!\cond
template <typename config_t>
concept deferred_config_element_concept = config_element_concept<config_t> &&
    std::is_base_of_v<deferred_config_element_base<std::remove_reference_t<config_t>>,
                                                   std::remove_reference_t<config_t>>;
//!\endcond

} // namespace seqan3::detail
