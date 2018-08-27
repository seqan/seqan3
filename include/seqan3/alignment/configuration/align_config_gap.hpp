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
 * \brief Provides gap configurations.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/utility.hpp>
#include <seqan3/alignment/gap/affine.hpp>
#include <seqan3/alignment/gap/linear.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3::detail
{
/*!\brief A configuration element for gaps.
 * \ingroup configuration
 *
 * \tparam gap_type The underlying gap class.
 */
template <typename gap_type>
struct align_config_gap
{
    //!\brief The actual value.
    gap_type value;
};

/*!\brief The gap adaptor enabling pipe notation.
 * \ingroup configuration
 *
 * \tparam gap_type A template-template class specifying the actual gap implementation.
 */
template <template <typename ...> typename gap_type>
//!\cond
    requires is_gap_config_v<gap_type<int>>
//!\endcond
struct align_config_gap_adaptor : public configuration_fn_base<align_config_gap_adaptor<gap_type>>
{

    /*!\brief Adds to the configuration a linear gap configuration element.
     * \param[in] cfg  The configuration to be extended.
     * \param[in] gs The gap score for the algorithm (usually negative).
     * \returns A new configuration containing the gap configuration element.
     */
    template <typename configuration_type,
              arithmetic_concept value_type>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_type>>
    //!\endcond
    constexpr auto invoke(configuration_type && cfg,
                          gap_score<value_type> const gs) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::gap, remove_cvref_t<configuration_type>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::gap));

        return std::forward<configuration_type>(cfg).push_front(align_config_gap<gap_type<value_type>>{{gs}});
    }

    /*!\brief Adds to the configuration an affine gap configuration element.
     * \param[in] cfg  The configuration to be extended.
     * \param[in] gs The gap score for the algorithm (usually negative).
     * \param[in] gos The gap open score for the algorithm (usually negative).
     * \returns A new configuration containing the gap configuration element.
     */
    template <typename configuration_type,
              arithmetic_concept value_type>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_type>>
    //!\endcond
    constexpr auto invoke(configuration_type && cfg,
                          gap_score<value_type> const gs,
                          gap_open_score<value_type> const gos) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::gap, remove_cvref_t<configuration_type>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::gap));

        return std::forward<configuration_type>(cfg).push_front(align_config_gap<gap_type<value_type>>{{gs, gos}});
    }
};

//!\brief Helper template meta-function associated with detail::align_config_gap.
//!\ingroup configuration
template <>
struct on_align_config<align_cfg::id::gap>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename is_type_specialisation_of<t, align_config_gap>::type;
};

//!\brief Mapping from the detail::align_config_gap type to its corresponding seqan3::align_cfg::id.
//!\ingroup configuration
template <typename value_type>
struct align_config_type_to_id<align_config_gap<value_type>>
{
    //!\brief The associated seqan3::align_cfg::id.
    static constexpr align_cfg::id value = align_cfg::id::gap;
};
} // namespace seqan3::detail

namespace seqan3::align_cfg
{
/*!\brief A configuration adaptor for linear gaps.
 * \ingroup configuration
 */
inline constexpr detail::align_config_gap_adaptor<seqan3::gap_linear> gap_linear;

/*!\brief A configuration adaptor for affine gaps.
 * \ingroup configuration
 */
inline constexpr detail::align_config_gap_adaptor<seqan3::gap_affine> gap_affine;

} // namespace seqan3::align_cfg
