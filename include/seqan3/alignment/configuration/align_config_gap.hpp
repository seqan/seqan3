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
 */

#pragma once

#include <seqan3/alignment/configuration/utility.hpp>
#include <seqan3/alignment/gap/linear.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3::detail
{
/*!\brief A configuration element for gaps.
 * \ingroup configuration
 *
 * \tparam gap_t The underlying gap class.
 */
template <typename gap_t>
class align_config_gap : public detail::config_element_base<align_config_gap<gap_t>>
{
    //!\brief Friend declaration to grant access to `value`.
    friend class detail::config_element_access<align_config_gap<gap_t>>;
    //!\brief The actual value.
    gap_t value;
};

/*!\brief The gap adaptor enabling pipe notation.
 * \ingroup configuration
 *
 * \tparam gap_t A template-template class specifying the actual gap implementation.
 */
template <template <typename ...> typename gap_t>
//!\cond
    requires is_gap_config_v<gap_t<int>>
//!\endcond
struct align_config_gap_adaptor : public configuration_fn_base<align_config_gap_adaptor<gap_t>>
{

    /*!\brief Adds to the configuration a gap configuration element.
     * \param[in] cfg  The configuration to be extended.
     * \param[in] cost The gap cost used to for the algorithm.
     * \returns A new configuration containing the gap configuration element.
     */
    template <typename configuration_t,
              typename value_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    constexpr auto invoke(configuration_t && cfg,
                          gap_cost<value_t> const cost) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::gap, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::gap));

        align_config_gap<gap_t<value_t>> tmp;
        tmp.get() = cost;
        return std::forward<configuration_t>(cfg).push_front(std::move(tmp));
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

//!\brief Mapping from the detail::align_config_gap type to it's corresponding seqan3::align_cfg::id.
//!\ingroup configuration
template <typename value_t>
struct align_config_type_to_id<align_config_gap<value_t>>
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

// inline constexpr detail::align_config_gap_adaptor<seqan3::gap_affine> gap_affine;
} // namespace seqan3::align_cfg
