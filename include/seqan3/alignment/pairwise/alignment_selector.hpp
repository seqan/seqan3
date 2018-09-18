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
 * \brief Provides seqan3::detail::alignment_selector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/type_list.hpp>

namespace seqan3::detail
{

/*!\brief Helper metafunction to determine the alignment result type based on the configuration.
 * \ingroup pairwise
 * \tparam seq1_t          The type of the first sequence.
 * \tparam seq2_t          The type of the second sequence.
 * \tparam configuration_t The configuration type. Must be of type seqan3::detail::configuration
 */
template <typename seq1_t, typename seq2_t, typename configuration_t>
struct determine_result_type
{
    //!\brief Helper function to determine the actual result type.
    static constexpr auto _determine()
    {
        using seq1_value_type = value_type_t<std::remove_reference_t<seq1_t>>;
        using seq2_value_type = value_type_t<std::remove_reference_t<seq2_t>>;
        using score_type      = int32_t;

        if constexpr (has_align_cfg_v<align_cfg::id::output, configuration_t>)
        {
            if constexpr (get<align_cfg::id::output>(configuration_t{}) ==
                               align_result_key::end)
                return align_result<type_list<uint32_t,
                                              score_type,
                                              std::pair<size_t, size_t>>>{};
            else if constexpr (get<align_cfg::id::output>(configuration_t{}) ==
                               align_result_key::begin)
                return align_result<type_list<uint32_t,
                                              score_type,
                                              std::pair<size_t, size_t>,
                                              std::pair<size_t, size_t>>>{};
            else if constexpr (get<align_cfg::id::output>(configuration_t{}) ==
                               align_result_key::trace)
                return align_result<type_list<uint32_t,
                                              score_type,
                                              std::pair<size_t, size_t>,
                                              std::pair<size_t, size_t>,
                                              std::tuple<std::vector<gapped<seq1_value_type>>,
                                                         std::vector<gapped<seq2_value_type>>>>>{};
            else
                return align_result<type_list<uint32_t,
                                              score_type>>{};
        }
        else
        {
            return align_result<type_list<uint32_t, score_type>>{};
        }
    }

    //!\brief The determined result type.
    using type = decltype(_determine());
};

/*!\brief Selects the correct alignment algorithm based on the algorithm configuration.
 * \ingroup pairwise
 * \tparam seq_tuple_t     A tuple like object containing the two source sequences.
 *                         Must model seqan3::tuple_like_concept.
 * \tparam configuration_t The specified configuration type.
 */
template <tuple_like_concept seq_tuple_t, typename configuration_t>
struct alignment_selector
{
    //!\brief The configuration stored globally for all alignment instances.
    configuration_t config;

    //!\brief The result type of invoking the algorithm.
    using result_type = typename determine_result_type<std::tuple_element_t<0, std::remove_reference_t<seq_tuple_t>>,
                                                       std::tuple_element_t<1, std::remove_reference_t<seq_tuple_t>>,
                                                       configuration_t>::type;

    /*!\brief Selects the corresponding alignment algorithm based on compile time and runtime decisions.
     * \param[in] seq The sequences as a tuple.
     * \returns A std::function object to be called within the alignment execution.
     *
     * \details
     *
     * \todo Write detail description and explain rationale about the function object.
     */
    template <tuple_like_concept _seq_tuple_t>
    auto select(_seq_tuple_t && seq)
    {
        //TODO Currently we only support edit_distance. We need would actually need real checks for this.
        std::function<result_type(result_type &)> func =
            pairwise_alignment_edit_distance_unbanded{std::get<0>(std::forward<_seq_tuple_t>(seq)),
                                                      std::get<1>(std::forward<_seq_tuple_t>(seq)),
                                                      config};
        return func;
    }
};
} // namespace seqan3::detail
