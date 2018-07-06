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
 * \brief Provides concepts for the configurator classes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>
#include <tuple>

#include <meta/meta.hpp>

#include <seqan3/core/algorithm/concept_pre.hpp>
#include <seqan3/std/concept/callable.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// configurator_concept
// ----------------------------------------------------------------------------
/*!\interface seqan3::detail::configurator_concept <>
 * \ingroup core_algorithm
 * \brief Concept for algorithm configurator.
 */

/*!\name Requirements for seqan3::detail::configurator_concept
 * \relates seqan3::detail::configurator_concept
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::detail::configurator_concept.
 * \{
 */
/*!\fn      get<I>(seqan3::configurator & cfg);
 * \brief   Returns the I-th configuration contained in seqan3::configurator.
 *
 * \tparam I   A non-type template argument to specify the position within the seqan3::configurator.
 * \param  cfg The seqan3::configurator to query.
 * \returns The stored configuration at the given position.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn      get<I>(seqan3::configurator && cfg);
 * \brief   Returns the I-th configuration contained in seqan3::configurator.
 *
 * \tparam I   A non-type template argument to specify the position within the seqan3::configurator.
 * \param  cfg The seqan3::configurator to query.
 * \returns The stored configuration at the given position.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn      get<I>(seqan3::configurator const & cfg);
 * \brief   Returns the I-th configuration contained in seqan3::configurator.
 *
 * \tparam I   A non-type template argument to specify the position within the seqan3::configurator.
 * \param  cfg The seqan3::configurator to query.
 * \returns The stored configuration at the given position.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn      get<I>(seqan3::configurator const && cfg);
 * \brief   Returns the I-th configuration contained in seqan3::configurator.
 *
 * \tparam I   A non-type template argument to specify the position within the seqan3::configurator.
 * \param  cfg The seqan3::configurator to query.
 * \returns The stored configuration at the given position.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn      get<type>(seqan3::configurator & cfg);
 * \brief   Returns the first configuration contained in seqan3::configurator which matches the given type.
 *
 * \tparam type The type for which the first match in seqan3::configurator should be returned.
 * \param  cfg  The seqan3::configurator to query.
 * \returns The stored configuration for the given type.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn      get<type>(seqan3::configurator && cfg);
 * \brief   Returns the first configuration contained in seqan3::configurator which matches the given type.
 *
 * \tparam type The type for which the first match in seqan3::configurator should be returned.
 * \param  cfg  The seqan3::configurator to query.
 * \returns The stored configuration for the given type.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn      get<type>(seqan3::configurator const & cfg);
 * \brief   Returns the first configuration contained in seqan3::configurator which matches the given type.
 *
 * \tparam type The type for which the first match in seqan3::configurator should be returned.
 * \param  cfg  The seqan3::configurator to query.
 * \returns The stored configuration for the given type.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn      get<type>(seqan3::configurator const && cfg);
 * \brief   Returns the first configuration contained in seqan3::configurator which matches the given type.
 *
 * \tparam type The type for which the first match in seqan3::configurator should be returned.
 * \param  cfg  The seqan3::configurator to query.
 * \returns The stored configuration for the given type.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\typedef typename seqan3::configurator::type_list_type type_list
 * \memberof seqan3::detail::configurator_concept
 * \brief Declares a seqan3::type_list over all contained configuration types.
 */
//!\}
//!\cond
template <typename configurator_t>
concept bool configurator_concept = requires (std::remove_reference_t<configurator_t> & cfg,
                                              std::remove_reference_t<configurator_t> const & c_cfg)
{
    // Expose the type_list.
    typename std::remove_reference_t<configurator_t>::type_list_type;

    // Support get functions.
    { get<0>(cfg) };
    { get<meta::at_c<typename std::remove_reference_t<configurator_t>::type_list_type, 0>>(cfg) };
    { get<0>(c_cfg) };
    { get<meta::at_c<typename std::remove_reference_t<configurator_t>::type_list_type, 0>>(c_cfg) };
    { get<0>(std::move(cfg)) };
    { get<meta::at_c<typename std::remove_reference_t<configurator_t>::type_list_type, 0>>(std::move(cfg)) };
    { get<0>(std::move(cfg)) };
    { get<meta::at_c<typename std::remove_reference_t<configurator_t>::type_list_type, 0>>(std::move(cfg)) };
};
//!\endcond

/*!\interface seqan3::detail::deferred_config_concept <>
 * \ingroup core_algorithm
 * \implements seqan3::detail::config_concept
 * \implements seqan3::invocable_concept
 * \brief Concept for an algorithm configuration.
 */

//!\cond
auto _callable_dummy = [] (auto && cfg)
{
    return std::forward<decltype(cfg)>(cfg);
};

template <typename config_t>
concept bool deferred_config_concept =
    config_concept<config_t> &&
    invocable_concept<std::remove_reference_t<config_t>,
                      decltype(_callable_dummy),
                      configurator<std::remove_reference_t<config_t>>>;
//!\endcond
} // namespace seqan3::detail
