// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::align_result_selector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <optional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_debug.hpp>
#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_traits/transformation_trait_or.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
/*!\brief Helper metafunction to select the alignment result type based on the configuration.
 * \ingroup pairwise_alignment
 * \tparam first_range_t   The type of the first sequence.
 * \tparam second_range_t  The type of the second sequence.
 * \tparam configuration_t The configuration type. Must be of type seqan3::detail::configuration
 */
template <std::ranges::forward_range first_range_t,
          std::ranges::forward_range second_range_t,
          typename configuration_t>
//!\cond
    requires is_type_specialisation_of_v<remove_cvref_t<configuration_t>, configuration>
//!\endcond
struct align_result_selector
{
private:
    //!\brief The user configured score type.
    using score_type = typename alignment_configuration_traits<
                                    std::remove_reference_t<configuration_t>>::original_score_type;

    //!\brief Helper function to determine the actual result type.
    static constexpr auto select()
    {
        static_assert(configuration_t::template exists<align_cfg::result>());

        if constexpr (configuration_t::template exists<align_cfg::result<with_back_coordinate_type, score_type>>())
        {
            return alignment_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate>{};
        }
        else if constexpr (configuration_t::template exists<align_cfg::result<with_front_coordinate_type, score_type>>())
        {
            return alignment_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate,
                                               alignment_coordinate>{};
        }
        else if constexpr (configuration_t::template exists<align_cfg::result<with_alignment_type, score_type>>())
        {
            // TODO: type `(first|second)_range_t &` is a hack to make the sequences viewable_ranges
            // https://github.com/seqan/seqan3/projects/10#card-33400557
            using alignment_type = typename make_pairwise_alignment_type<first_range_t &, second_range_t &>::type;

            return alignment_result_value_type<uint32_t,
                                               score_type,
                                               alignment_coordinate,
                                               alignment_coordinate,
                                               alignment_type>{};
        }
        else
        {
            return alignment_result_value_type<uint32_t, score_type>{};
        }
    }

    //!\brief Augments the result type with the debug output for score and trace matrix.
    template <typename alignment_result_value_t>
    static constexpr auto augment_if_debug(alignment_result_value_t)
    {
        if constexpr (configuration_t::template exists<detail::debug_mode>())
        {
            using as_type_list = transfer_template_args_onto_t<alignment_result_value_t, type_list>;
            using score_matrix_t = two_dimensional_matrix<std::optional<score_type>,
                                                          std::allocator<std::optional<score_type>>,
                                                          matrix_major_order::column>;
            using trace_matrix_t = two_dimensional_matrix<std::optional<trace_directions>,
                                                          std::allocator<std::optional<trace_directions>>,
                                                          matrix_major_order::column>;

            if constexpr (configuration_t::template exists<align_cfg::result<with_alignment_type>>())
            {
                using with_score_t = list_traits::replace_at<score_matrix_t, 5, as_type_list>;
                return transfer_template_args_onto_t<list_traits::replace_at<trace_matrix_t, 6, with_score_t>,
                                                     alignment_result_value_type>{};
            }
            else
            {
                return transfer_template_args_onto_t<list_traits::replace_at<score_matrix_t, 5, as_type_list>,
                                                     alignment_result_value_type>{};
            }
        }
        else  // Return as is.
        {
            return alignment_result_value_t{};
        }
    }

public:
    //!\brief The selected result type.
    using type = decltype(augment_if_debug(select()));
};

} // namespace seqan3::detail
