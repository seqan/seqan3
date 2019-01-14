// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/range/view/persist.hpp>

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
        using seq1_value_type = gapped<value_type_t<std::remove_reference_t<seq1_t>>>;
        using seq2_value_type = gapped<value_type_t<std::remove_reference_t<seq2_t>>>;
        using score_type      = int32_t;

        if constexpr (std::remove_reference_t<configuration_t>::template exists<align_cfg::result>())
        {
            if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                    with_end_position_type>)
                return align_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate>{};
            else if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                         with_begin_position_type>)
                return align_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate,
                                               alignment_coordinate>{};
            else if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                         with_trace_type>)
                return align_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate,
                                               alignment_coordinate,
                                               std::pair<std::vector<seq1_value_type>,
                                                         std::vector<seq2_value_type>>>{};
            else
                return align_result_value_type<uint32_t, score_type>{};
        }
        else
        {
            return align_result_value_type<uint32_t, score_type>{};
        }
    }

    //!\brief The determined result type.
    using type = align_result<decltype(_determine())>;
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
            pairwise_alignment_edit_distance_unbanded{std::get<0>(std::forward<_seq_tuple_t>(seq)) | view::persist,
                                                      std::get<1>(std::forward<_seq_tuple_t>(seq)) | view::persist,
                                                      config};
        return func;
    }
};
} // namespace seqan3::detail
