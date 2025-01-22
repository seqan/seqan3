// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides helper type traits for the configuration and execution of the alignment algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_debug.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_on_result.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/alignment/configuration/align_config_result_type.hpp>
#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_vectorised.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/detail/bits_of.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/type_traits/function_traits.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

//------------------------------------------------------------------------------
// chunked_indexed_sequence_pairs
//------------------------------------------------------------------------------

/*!\brief A transformation trait to retrieve the chunked range over indexed sequence pairs.
 * \ingroup alignment_pairwise
 * \implements seqan3::transformation_trait
 *
 * \tparam sequence_pairs_t The type of the sequences to be transformed; must model seqan3::detail::sequence_pair_range.
 *
 * \details
 *
 * This transformation trait transforms a range over sequence pairs into a range over indexed sequence pairs.
 * In addition, the range is chunked which is the common interface for alignment algorithms.
 * The returned type models seqan3::detail::indexed_sequence_pair_range.
 */
template <typename sequence_pairs_t>
    requires sequence_pair_range<std::remove_reference_t<sequence_pairs_t>>
struct chunked_indexed_sequence_pairs
{
    //!\brief The transformed type that models seqan3::detail::indexed_sequence_pair_range.
    using type = decltype(views::zip(std::declval<sequence_pairs_t>(), std::views::iota(0)) | views::chunk(1));
};

//------------------------------------------------------------------------------
// alignment_configuration_traits
//------------------------------------------------------------------------------

/*!\brief A traits type for the alignment algorithm that exposes static information stored within the alignment
 *        configuration object.
 * \ingroup alignment_pairwise
 *
 * \tparam configuration_t The type of the alignment configuration object; must be a specialisation of
 *         seqan3::configuration.
 */
template <typename configuration_t>
    requires is_type_specialisation_of_v<std::remove_cv_t<configuration_t>, configuration>
struct alignment_configuration_traits
{
private:
    /*!\brief An index type (i.e. unsigned integral) for a score_type which has the same bit size.
     * \tparam score_t The score type for which the corresponding index type shall be returned; must model
     *                 seqan3::arithmetic.
     *
     * \details
     *
     * If the score type models std::integral it will simply be converted to its unsigned counterpart.
     * For floating point types the bit size is determined and the corresponding minimal viable
     * unsigned integral type (seqan3::detail::min_viable_uint_t) is returned.
     */
    template <arithmetic score_t>
    using select_scalar_index_t = min_viable_uint_t<1ull << (bits_of<score_t> - 1)>;

    //!\brief Helper function to determine the alignment result type.
    static constexpr auto determine_alignment_result_type() noexcept
    {
        if constexpr (configuration_t::template exists<align_cfg::detail::result_type>())
        {
            using result_type_cfg_t = std::remove_cvref_t<decltype(seqan3::get<align_cfg::detail::result_type>(
                std::declval<configuration_t>()))>;
            return typename result_type_cfg_t::type{}; // Access the stored result_type.
        }
        else
        {
            return empty_type{};
        }
    }

public:
    //!\brief Flag to indicate vectorised mode.
    static constexpr bool is_vectorised = configuration_t::template exists<align_cfg::vectorised>();
    //!\brief Flag indicating whether parallel alignment mode is enabled.
    static constexpr bool is_parallel = configuration_t::template exists<align_cfg::parallel>();
    //!\brief Flag indicating whether global alignment method is enabled.
    static constexpr bool is_global = configuration_t::template exists<seqan3::align_cfg::method_global>();
    //!\brief Flag indicating whether local alignment mode is enabled.
    static constexpr bool is_local = configuration_t::template exists<seqan3::align_cfg::method_local>();
    //!\brief Flag indicating whether banded alignment mode is enabled.
    static constexpr bool is_banded = configuration_t::template exists<align_cfg::band_fixed_size>();
    //!\brief Flag indicating whether debug mode is enabled.
    static constexpr bool is_debug = configuration_t::template exists<detail::debug_mode>();
    //!\brief Flag indicating whether a user provided callback was given.
    static constexpr bool is_one_way_execution = configuration_t::template exists<align_cfg::on_result>();
    //!\brief The selected scoring scheme.
    using scoring_scheme_type = decltype(get<align_cfg::scoring_scheme>(std::declval<configuration_t>()).scheme);
    //!\brief The alphabet of the selected scoring scheme.
    using scoring_scheme_alphabet_type = typename scoring_scheme_type::alphabet_type;
    //!\brief The original score type selected by the user.
    using original_score_type = typename std::remove_reference_t<decltype(std::declval<configuration_t>().get_or(
        align_cfg::score_type<int32_t>{}))>::type;
    //!\brief The score type for the alignment algorithm.
    using score_type = std::conditional_t<is_vectorised, simd_type_t<original_score_type>, original_score_type>;
    //!\brief The trace directions type for the alignment algorithm.
    using trace_type = std::conditional_t<is_vectorised, simd_type_t<original_score_type>, trace_directions>;
    //!\brief The alignment result type if present. Otherwise seqan3::detail::empty_type.
    using alignment_result_type = decltype(determine_alignment_result_type());
    //!\brief The type of the matrix index.
    using matrix_index_type =
        std::conditional_t<is_vectorised, simd_type_t<select_scalar_index_t<original_score_type>>, size_t>;
    //!\brief The type of the matrix coordinate.
    using matrix_coordinate_type =
        lazy_conditional_t<is_vectorised, lazy<simd_matrix_coordinate, matrix_index_type>, matrix_coordinate>;

    //!\brief The number of alignments that can be computed in one simd vector.
    static constexpr size_t alignments_per_vector = []() constexpr
    {
        if constexpr (is_vectorised)
            return simd_traits<score_type>::length;
        else
            return 1;
    }();
    //!\brief Flag indicating whether the score shall be computed.
    static constexpr bool compute_score = configuration_t::template exists<align_cfg::output_score>();
    //!\brief Flag indicating whether the end positions shall be computed.
    static constexpr bool compute_end_positions = configuration_t::template exists<align_cfg::output_end_position>();
    //!\brief Flag indicating whether the begin positions shall be computed.
    static constexpr bool compute_begin_positions =
        configuration_t::template exists<align_cfg::output_begin_position>();
    //!\brief Flag indicating whether the sequence alignment shall be computed.
    static constexpr bool compute_sequence_alignment = configuration_t::template exists<align_cfg::output_alignment>();
    //!\brief Flag indicating whether the id of the first sequence shall be returned.
    static constexpr bool output_sequence1_id = configuration_t::template exists<align_cfg::output_sequence1_id>();
    //!\brief Flag indicating whether the id of the second sequence shall be returned.
    static constexpr bool output_sequence2_id = configuration_t::template exists<align_cfg::output_sequence2_id>();
    //!\brief Flag indicating if any output option was set.
    static constexpr bool has_output_configuration = compute_score || compute_end_positions || compute_begin_positions
                                                  || compute_sequence_alignment || output_sequence1_id
                                                  || output_sequence2_id;
    //!\brief Flag indicating whether the trace matrix needs to be computed.
    static constexpr bool requires_trace_information = compute_begin_positions || compute_sequence_alignment;
};

//------------------------------------------------------------------------------
// alignment_function_traits
//------------------------------------------------------------------------------

/*!\brief A traits class to provide a uniform access to the properties of the wrapped alignment algorithm.
 * \ingroup alignment_pairwise
 *
 * \tparam function_t The type of the std::function object that stores the alignment algorithm as the target.
 */
template <typename function_t>
struct alignment_function_traits
{
    //!\brief The type of the sequence input to the alignment algorithm.
    using sequence_input_type = typename function_traits<function_t>::template argument_type_at<0>;
    //!\brief The type of the callback function called when a result was computed.
    using callback_type = typename function_traits<function_t>::template argument_type_at<1>;
    //!\brief The type of the alignment result to be computed.
    using alignment_result_type = typename function_traits<callback_type>::template argument_type_at<0>;
};

} // namespace seqan3::detail
