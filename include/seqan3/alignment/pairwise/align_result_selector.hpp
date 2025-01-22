// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::align_result_selector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <optional>
#include <ranges>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_debug.hpp>
#include <seqan3/alignment/decorator/gap_decorator.hpp>
#include <seqan3/alignment/matrix/detail/advanceable_alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/type_traits/basic.hpp>
#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{
/*!\brief Helper metafunction to select the alignment result type based on the configuration.
 * \ingroup alignment_pairwise
 * \tparam first_range_t   The type of the first sequence.
 * \tparam second_range_t  The type of the second sequence.
 * \tparam configuration_t The configuration type. Must be of type seqan3::detail::configuration
 */
template <std::ranges::forward_range first_range_t, std::ranges::forward_range second_range_t, typename configuration_t>
    requires is_type_specialisation_of_v<std::remove_cvref_t<configuration_t>, configuration>
struct align_result_selector
{
private:
    //!\brief The traits type used for the alignment.
    using traits_type = alignment_configuration_traits<configuration_t>;
    //!\brief The user configured score type.
    using score_type =
        typename alignment_configuration_traits<std::remove_reference_t<configuration_t>>::original_score_type;
    //!\brief The type to indicate that an output option was not configured.
    using disabled_type = std::nullopt_t *;
    //!\brief Score matrix type in debug mode.
    using debug_score_matrix_type = two_dimensional_matrix<std::optional<score_type>,
                                                           std::allocator<std::optional<score_type>>,
                                                           matrix_major_order::column>;
    //!\brief Trace matrix type in debug mode.
    using debug_trace_matrix_type = two_dimensional_matrix<std::optional<trace_directions>,
                                                           std::allocator<std::optional<trace_directions>>,
                                                           matrix_major_order::column>;
    //!\brief The configured score type if selected.
    using configured_score_type = std::conditional_t<traits_type::compute_score, score_type, disabled_type>;
    //!\brief The configured end position type if selected.
    using configured_end_position_type = std::conditional_t<traits_type::compute_end_positions,
                                                            seqan3::detail::advanceable_alignment_coordinate<>,
                                                            disabled_type>;
    //!\brief The configured begin position type if selected.
    using configured_begin_position_type = std::conditional_t<traits_type::compute_begin_positions,
                                                              seqan3::detail::advanceable_alignment_coordinate<>,
                                                              disabled_type>;
    //!\brief The configured alignment type if selected.
    using configured_alignment_type =
        typename lazy_conditional_t<traits_type::compute_sequence_alignment,
                                    lazy<make_pairwise_alignment_type, first_range_t &, second_range_t &>,
                                    std::type_identity<disabled_type>>::type;

    //!\brief The configured sequence id type for the first sequence if selected.
    using configured_sequence1_id_type = std::conditional_t<traits_type::output_sequence1_id, uint32_t, disabled_type>;
    //!\brief The configured sequence id type for the second sequence if selected.
    using configured_sequence2_id_type = std::conditional_t<traits_type::output_sequence2_id, uint32_t, disabled_type>;

    //!\brief The debug score matrix type if selected.
    using configured_debug_score_matrix_type =
        std::conditional_t<traits_type::is_debug, debug_score_matrix_type, disabled_type>;

    //!\brief The debug trace matrix type if selected.
    using configured_debug_trace_matrix_type =
        std::conditional_t<traits_type::is_debug && traits_type::compute_sequence_alignment,
                           debug_trace_matrix_type,
                           disabled_type>;

public:
    //!\brief The selected result type.
    using type = alignment_result_value_type<configured_sequence1_id_type,
                                             configured_sequence2_id_type,
                                             configured_score_type,
                                             configured_end_position_type,
                                             configured_begin_position_type,
                                             configured_alignment_type,
                                             configured_debug_score_matrix_type,
                                             configured_debug_trace_matrix_type>;
};

} // namespace seqan3::detail
