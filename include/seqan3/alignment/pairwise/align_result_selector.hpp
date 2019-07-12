// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
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
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/type_traits/transformation_trait_or.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A helper class to define the alignment return type.
 * \tparam first_t  Type of the first sequence.
 * \tparam second_t Type of the second sequence.
 * \details
 * The type uses the gap decorator if RandomAccessRange and SizedRange are met for both input sequences.
 */
template <typename first_t, typename second_t>
struct alignment_type;

//!\overload
template <typename first_t, typename second_t>
    requires std::ranges::RandomAccessRange<first_t> &&
             std::ranges::SizedRange<first_t> &&
             std::ranges::RandomAccessRange<second_t> &&
             std::ranges::SizedRange<second_t>
struct alignment_type<first_t, second_t>
{
    //!\brief The alignment type with gap decorator.
    using type = std::tuple<gap_decorator<all_view<first_t &>>, gap_decorator<all_view<second_t &>>>;
};

/*!\brief Helper metafunction to select the alignment result type based on the configuration.
 * \ingroup pairwise_alignment
 * \tparam first_bach_t    The type of the first sequence.
 * \tparam second_range_t  The type of the second sequence.
 * \tparam configuration_t The configuration type. Must be of type seqan3::detail::configuration
 */
template <std::ranges::ForwardRange first_range_t,
          std::ranges::ForwardRange second_range_t,
          typename configuration_t>
//!\cond
    requires is_type_specialisation_of_v<remove_cvref_t<configuration_t>, configuration>
//!\endcond
struct align_result_selector
{
    //!\brief Helper function to determine the actual result type.
    static constexpr auto _determine()
    {
        using score_type            = int32_t;

        if constexpr (std::remove_reference_t<configuration_t>::template exists<align_cfg::result>())
        {
            if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                    with_back_coordinate_type>)
            {
                return alignment_result_value_type<uint32_t,
                                                   score_type,
                                                   alignment_coordinate>{};
            }
            else if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                         with_front_coordinate_type>)
            {
                return alignment_result_value_type<uint32_t,
                                                   score_type,
                                                   alignment_coordinate,
                                                   alignment_coordinate>{};
            }
            else if constexpr (std::Same<remove_cvref_t<decltype(get<align_cfg::result>(configuration_t{}).value)>,
                                         with_alignment_type>)
            {
                // Due to an error with gcc8 we define these types beforehand.
                using first_gapped_seq_type = gapped<value_type_t<first_range_t>>;
                using second_gapped_seq_type = gapped<value_type_t<second_range_t>>;

                // We use vectors of gapped sequence if the gap decorator cannot be used.
                using fallback_t = std::tuple<std::vector<first_gapped_seq_type>, std::vector<second_gapped_seq_type>>;

                // If the ranges are RandomAccess and Sized we can use the Gap Decorator, otherwise fallback_t.
                using decorator_t = alignment_type<first_range_t, second_range_t>;

                return alignment_result_value_type<uint32_t,
                                                   score_type,
                                                   alignment_coordinate,
                                                   alignment_coordinate,
                                                   detail::transformation_trait_or_t<decorator_t, fallback_t>>{};
            }
            else
            {
                return alignment_result_value_type<uint32_t, score_type>{};
            }
        }
        else
        {
            return alignment_result_value_type<uint32_t, score_type>{};
        }
    }

    //!\brief The determined result type.
    using type = decltype(_determine());
};

} // namespace seqan3::detail
