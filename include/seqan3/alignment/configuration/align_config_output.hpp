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
 * \brief Provides configuration for alignment output.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/utility.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3::detail
{
/*!\brief A configuration element for global alignment.
 * \ingroup configuration
 */
template <align_result_key e>
struct align_config_output
{
    //!\brief The value of align_config_output.
    align_result_key value{e};
};

/*!\brief The output adaptor enabling pipe notation.
 * \ingroup configuration
 */
template <align_result_key e>
struct align_config_output_adaptor : public configuration_fn_base<align_config_output_adaptor<e>>
{
    static_assert(e != align_result_key::id,
                  "The id field is only usable to recover the id of the alignment within the alignment result. "
                  "You need to at least use align_result_key::score for the alignment configuration.");
    /*!\brief Adds to the configuration an output configuration element.
     * \param[in] cfg  The configuration to be extended.
     * \tparam configuration_t The type of the underlying configuration scheme.
     *                         Is required to fulfill the seqan3::detail::is_algorithm_configuration requirement.
     * \returns A new configuration containing the output configuration element.
     */
    template <typename configuration_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    constexpr auto invoke(configuration_t && cfg) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::output, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::output));

        return std::forward<configuration_t>(cfg).push_front(align_config_output<e>{});
    }
};

//!\brief Helper template meta-function associated with seqan3::detail::align_config_output.
//!\ingroup configuration
template <>
struct on_align_config<align_cfg::id::output>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename is_value_specialisation_of<t, align_config_output>::type;
};

//!\brief Mapping from the seqan3::detail::align_config_output type to its corresponding seqan3::align_cfg::id.
//!\ingroup configuration
template <align_result_key e>
struct align_config_type_to_id<align_config_output<e>>
{
    //!\brief The associated seqan3::align_cfg::id.
    static constexpr align_cfg::id value = align_cfg::id::output;
};
} // namespace seqan3::detail

namespace seqan3::align_cfg
{
/*!\brief A configuration adaptor for alignment output.
 * \ingroup configuration
 */
template <align_result_key e>
inline constexpr detail::align_config_output_adaptor<e> output;
} // namespace seqan3::align_cfg
