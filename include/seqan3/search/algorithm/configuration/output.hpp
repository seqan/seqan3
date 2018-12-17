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
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/search/algorithm/configuration/detail.hpp>

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

/*!\brief Configuration element to determine the output type of hits.
 * \ingroup search_configuration
 */
template <typename output_t>
//!\cond
    requires std::Same<remove_cvref_t<output_t>, detail::search_output_text_position> ||
             std::Same<remove_cvref_t<output_t>, detail::search_output_index_iterator>
//!\endcond
struct output : public pipeable_config_element<output<output_t>, output_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output};
};

/*!\name Type deduction guides
 * \relates seqan3::search_cfg::output
 * \{
 */

//!\brief Deduces search output type from constructor argument.
template <typename output_t>
output(output_t) -> output<remove_cvref_t<output_t>>;
//!\}

} // namespace seqan3::search_cfg

//!\}
