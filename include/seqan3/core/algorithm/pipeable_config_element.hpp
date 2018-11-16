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
 * \brief Provides seqan3::pipeable_config_element.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/metafunction/basic.hpp>

namespace seqan3
{

/*!\brief Adds pipe interface to configuration elements.
 * \ingroup algorithm
 */
struct pipeable_config_element
{};

/*!\name Pipe operator
 * \{
 */

/*!\brief Combines left operand with right operand in a new seqan3::configuration.
 * \relates seqan3::pipeable_config_element
 * \ingroup algorithm
 *
 * \param[in] lhs The left operand. Must satisfy seqan3::detail::config_element_concept.
 * \param[in] rhs The right operand. Must satisfy seqan3::detail::config_element_concept.
 *
 * \returns The new configuration instance with containing both configuration elements.
 *
 * \details
 *
 * Effectively calls seqan3::configuration::push_back on a temporary created seqan3::configuration from the
 * left operand.
 */
template <typename lhs_t, typename rhs_t>
//!\cond
    requires std::is_base_of_v<pipeable_config_element, remove_cvref_t<lhs_t>> &&
             std::is_base_of_v<pipeable_config_element, remove_cvref_t<rhs_t>>
//!\endcond
constexpr auto operator|(lhs_t && lhs, rhs_t && rhs)
{
    return configuration{std::forward<lhs_t>(lhs)}.push_back(std::forward<rhs_t>(rhs));
}

/*!\brief Combines left operand with right operand in a new seqan3::configuration.
 * \relates seqan3::pipeable_config_element
 * \ingroup algorithm
 *
 * \param[in] lhs The left operand. Must be of type seqan3::configuration.
 * \param[in] rhs The right operand. Must satisfy seqan3::detail::config_element_concept.
 *
 * \returns The new configuration instance with containing the previous configuration elements and the new one.
 *
 * \details
 *
 * Effectively calls seqan3::configuration::push_back on the left operand.
 */
template <typename lhs_t, typename rhs_t>
//!\cond
    requires detail::is_algorithm_configuration_v<remove_cvref_t<lhs_t>> &&
             std::is_base_of_v<pipeable_config_element, remove_cvref_t<rhs_t>>
//!\endcond
constexpr auto operator|(lhs_t && lhs, rhs_t && rhs)
{
    return std::forward<lhs_t>(lhs).push_back(std::forward<rhs_t>(rhs));
}
//!\}
} // namespace seqan3
