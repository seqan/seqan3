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
template <typename derived_t, typename value_type>
struct pipeable_config_element
{
    //!\brief The stored config value.
    value_type value;

    /*!\name Pipe operator
     * \{
     */
    /*!\brief Combines two configuration elements to a seqan3::configuration.
     * \tparam other_derived_t The derived type of the right hand side operand.
     * \tparam other_value_t   The value type of the right hand side operand.
     * \param[in] lhs          The left hand operand.
     * \param[in] rhs          The right hand operand.
     * \returns A new seqan3::configuration containing `lhs` and `rhs`.
     *
     * \details
     *
     * Effectively calls seqan3::configuration::push_back on a temporary created seqan3::configuration from the
     * left operand.
     */
    template <typename other_derived_t, typename other_value_t>
    friend constexpr auto operator|(pipeable_config_element && lhs,
                                    pipeable_config_element<other_derived_t, other_value_t> && rhs)
    {
        return configuration{static_cast<derived_t &&>(lhs)}.push_back(static_cast<other_derived_t &&>(rhs));
    }

    //!\overload
    template <typename other_derived_t, typename other_value_t>
    friend constexpr auto operator|(pipeable_config_element && lhs,
                                    pipeable_config_element<other_derived_t, other_value_t> const & rhs)
    {
        return configuration{static_cast<derived_t &&>(lhs)}.push_back(static_cast<other_derived_t const &>(rhs));
    }

    //!\overload
    template <typename other_derived_t, typename other_value_t>
    friend constexpr auto operator|(pipeable_config_element const & lhs,
                                    pipeable_config_element<other_derived_t, other_value_t> && rhs)
    {
        return configuration{static_cast<derived_t const &>(lhs)}.push_back(static_cast<other_derived_t &&>(rhs));
    }

    //!\overload
    template <typename other_derived_t, typename other_value_t>
    friend constexpr auto operator|(pipeable_config_element const & lhs,
                                    pipeable_config_element<other_derived_t, other_value_t> const & rhs)
    {
        return configuration{static_cast<derived_t const &>(lhs)}.push_back(static_cast<other_derived_t const &>(rhs));
    }

    /*!\brief Combines two configuration elements to a seqan3::configuration.
     * \tparam configs_t  A template parameter pack for the given seqan3::configuration.
     * \param[in] lhs     The left hand operand.
     * \param[in] rhs     The right hand operand.
     * \returns A new seqan3::configuration adding `rhs` to the passed `lhs` object.
     *
     * \details
     *
     * Effectively calls seqan3::configuration::push_back on the left operand.
     */
    template <typename ...configs_t>
    friend constexpr auto operator|(configuration<configs_t...> && lhs,
                                    pipeable_config_element && rhs)
    {
        return std::move(lhs).push_back(static_cast<derived_t &&>(rhs));
    }

    //!\overload
    template <typename ...configs_t>
    friend constexpr auto operator|(configuration<configs_t...> const & lhs,
                                    pipeable_config_element && rhs)
    {
        return lhs.push_back(static_cast<derived_t &&>(rhs));
    }

    //!\overload
    template <typename ...configs_t>
    friend constexpr auto operator|(configuration<configs_t...> && lhs,
                                    pipeable_config_element const & rhs)
    {
        return std::move(lhs).push_back(static_cast<derived_t const &>(rhs));
    }

    //!\overload
    template <typename ...configs_t>
    friend constexpr auto operator|(configuration<configs_t...> const & lhs,
                                    pipeable_config_element const & rhs)
    {
        return lhs.push_back(static_cast<derived_t const &>(rhs));
    }
    //!\}
};

} // namespace seqan3
