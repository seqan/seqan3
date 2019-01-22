// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::align_result_selector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Helper metafunction to select the alignment result type based on the configuration.
 * \ingroup pairwise_alignment
 * \tparam first_bach_t    The type of the first sequence.
 * \tparam second_batch_t  The type of the second sequence.
 * \tparam configuration_t The configuration type. Must be of type seqan3::detail::configuration
 */
template <std::ranges::ForwardRange first_batch_t,
          std::ranges::ForwardRange second_batch_t,
          typename configuration_t>
//!\cond
    requires is_type_specialisation_of_v<remove_cvref_t<configuration_t>, configuration>
//!\endcond
struct align_result_selector
{
    //!\brief Helper function to determine the actual result type.
    static constexpr auto _determine()
    {
        using first_seq_value_type  = gapped<value_type_t<first_batch_t>>;
        using second_seq_value_type = gapped<value_type_t<second_batch_t>>;
        using score_type            = int32_t;

        if constexpr (std::remove_reference_t<configuration_t>::template exists<align_cfg::result>())
        {
            if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                    with_end_position_type>)
            {
                return align_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate>{};
            }
            else if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                         with_begin_position_type>)
            {
                return align_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate,
                                               alignment_coordinate>{};
            }
            else if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                         with_trace_type>)
            {
                return align_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate,
                                               alignment_coordinate,
                                               std::tuple<std::vector<first_seq_value_type>,
                                                          std::vector<second_seq_value_type>>>{};
            }
            else
            {
                return align_result_value_type<uint32_t, score_type>{};
            }
        }
        else
        {
            return align_result_value_type<uint32_t, score_type>{};
        }
    }

    //!\brief The determined result type.
    using type = decltype(_determine());
};
} // namespace seqan3::detail
