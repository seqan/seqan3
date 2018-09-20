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
 * \brief Provides the configuration for returning positions in the text.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/search/algorithm/configuration/utility.hpp>

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

//!\brief Type for the "index_iterator" value for the configuration element "output".
struct search_output_index_iterator {};
//!\brief Type for the "text_position" value for the configuration element "output".
struct search_output_text_position {};

} // namespace seqan3::detail

namespace seqan3::search_cfg
{

//!\brief Configuration element to receive all hits within the error bounds.
inline detail::search_output_index_iterator constexpr index_iterator;
//!\brief Configuration element to receive all hits within the lowest number of errors.
inline detail::search_output_text_position constexpr text_position;

} // namespace seqan3::search_cfg

namespace seqan3::detail
{
/*!\brief Configuration element to determine the output type of hits.
 * \ingroup search_configuration
 */
template <typename output_t>
struct search_config_output
{
    //!\cond
    output_t value;
    //!\endcond
};

/*!\brief The seqan3::search_cfg::output adaptor enabling pipe notation.
 * \ingroup search_configuration
 */
template <template <typename ...> typename search_config_output_type>
struct search_config_output_adaptor :
    public configuration_fn_base<search_config_output_adaptor<search_config_output_type>>
{

    /*!\brief Adds to the configuration the seqan3::search_cfg::output configuration element.
     * \param[in] cfg The configuration to be extended.
     * \param[in] output The output type determining the type of the reported hits.
     * \returns A new configuration containing the seqan3::search_cfg::output configuration element.
     */
    template <typename configuration_t, typename output_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>> &&
                 (std::Same<remove_cvref_t<output_t>, search_output_index_iterator> ||
                  std::Same<remove_cvref_t<output_t>, search_output_text_position>)
    //!\endcond
    constexpr auto invoke(configuration_t && cfg, output_t output) const
    {
        static_assert(is_valid_search_configuration_v<search_cfg::id::output, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(search_cfg::id::output));

        return std::forward<configuration_t>(cfg).push_front(search_config_output<output_t>{std::move(output)});
    }
};

//!\brief Helper template meta-function associated with detail::search_config_output.
//!\ingroup search_configuration
template <>
struct on_search_config<search_cfg::id::output>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename is_type_specialisation_of<t, search_config_output>::type;
};

//!\brief Mapping from the detail::search_config_output type to it's corresponding seqan3::search_cfg::id.
//!\ingroup search_configuration
template <typename output_t>
struct search_config_type_to_id<search_config_output<output_t>>
{
    //!\brief The associated seqan3::search_cfg::id.
    static constexpr search_cfg::id value = search_cfg::id::output;
};
} // namespace seqan3::detail

namespace seqan3::search_cfg
{
/*!\brief Configuration element to determine the output type of hits.
 * \ingroup search_configuration
 */
inline detail::search_config_output_adaptor<seqan3::detail::search_config_output> constexpr output;

} // namespace seqan3::search_cfg

//!\}
