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
 * \brief Provides scoring configurations.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/utility.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3::detail
{
/*!\brief A configuration element for alignment scoring.
 * \ingroup configuration
 * \tparam scoring_scheme_type The type of the underlying scoring scheme.
 */
template <typename scoring_scheme_type>
struct align_config_score
{
    //!\brief Holds the actual scoring scheme.
    scoring_scheme_type value;
};

/*!\brief The score adaptor enabling pipe notation.
 * \ingroup configuration
 */
struct align_config_score_adaptor : public configuration_fn_base<align_config_score_adaptor>
{
    /*!\brief Adds to the configuration a score configuration element.
     * \tparam configuration_type The type of the configuration to be extended.
     * \tparam scoring_scheme_type The type of the scoring scheme that is included in the new configuration.
     * \param[in] cfg  The configuration to be extended.
     * \param[in] scheme The scoring scheme for the algorithm.
     * \returns A new configuration containing the score configuration element.
     */
    template <typename configuration_type, typename scoring_scheme_type>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_type>>
    //!\endcond
    constexpr auto invoke(configuration_type && cfg, scoring_scheme_type const scheme) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::score, remove_cvref_t<configuration_type>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::score));

        return std::forward<configuration_type>(cfg).push_front(align_config_score<scoring_scheme_type>
                                                                {std::move(scheme)});
    }
};

/*!\brief Helper template meta-function associated with seqan3::detail::align_config_score.
 * \ingroup configuration
 */
template <>
struct on_align_config<align_cfg::id::score>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename is_type_specialisation_of<t, align_config_score>::type;
};

/*!\brief Mapping from the seqan3::detail::align_config_score type to its corresponding seqan3::align_cfg::id.
 * \tparam scoring_scheme_type The type of the scoring scheme in the configuration element.
 * \ingroup configuration
 */
template <typename scoring_scheme_type>
struct align_config_type_to_id<align_config_score<scoring_scheme_type>>
{
    //!\brief The associated seqan3::align_cfg::id.
    static constexpr align_cfg::id value = align_cfg::id::score;
};

} // namespace seqan3::detail

namespace seqan3::align_cfg
{

/*!\brief A configuration adaptor for alignment scoring.
 * \ingroup configuration
 */
inline constexpr detail::align_config_score_adaptor score;

} // namespace seqan3::align_cfg
