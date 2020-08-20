// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides helper type traits for the configuration and execution of the alignment algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_result_type.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_debug.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_on_result.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_vectorised.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/type_traits/function.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

//------------------------------------------------------------------------------
// chunked_indexed_sequence_pairs
//------------------------------------------------------------------------------

/*!\brief A transformation trait to retrieve the chunked range over indexed sequence pairs.
 * \ingroup pairwise_alignment
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
//!\cond
    requires sequence_pair_range<std::remove_reference_t<sequence_pairs_t>>
//!\endcond
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
 * \ingroup pairwise_alignment
 *
 * \tparam configuration_t The type of the alignment configuration object; must be a specialisation of
 *         seqan3::configuration.
 */
template <typename configuration_t>
//!\cond
    requires is_type_specialisation_of_v<std::remove_cv_t<configuration_t>, configuration>
//!\endcond
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
    using select_scalar_index_t = min_viable_uint_t<1ull << (sizeof_bits<score_t> - 1)>;

    //!\brief Helper function to determine the alignment result type.
    static constexpr auto determine_alignment_result_type() noexcept
    {
        if constexpr (configuration_t::template exists<result_type_tag>())
        {
            using wrapped_result_t =
                decltype(seqan3::get<result_type_tag>(std::declval<configuration_t>()).value);
            return typename wrapped_result_t::type{};  // Unwrap the type_identity.
        }
        else
        {
            return empty_type{};
        }
    }

public:
    //!\brief Flag to indicate vectorised mode.
    static constexpr bool is_vectorised =
        configuration_t::template exists<remove_cvref_t<decltype(align_cfg::vectorised)>>();
    //!\brief Flag indicating whether parallel alignment mode is enabled.
    static constexpr bool is_parallel = configuration_t::template exists<align_cfg::parallel>();
    //!\brief Flag indicating whether global alignment method is enabled.
    static constexpr bool is_global =
        configuration_t::template exists<seqan3::align_cfg::method_global>();
    //!\brief Flag indicating whether global alignment mode with free ends is enabled.
    static constexpr bool with_free_end_gaps = configuration_t::template exists<align_cfg::aligned_ends>();
    //!\brief Flag indicating whether local alignment mode is enabled.
    static constexpr bool is_local = configuration_t::template exists<seqan3::detail::method_local_tag>();
    //!\brief Flag indicating whether banded alignment mode is enabled.
    static constexpr bool is_banded = configuration_t::template exists<align_cfg::band_fixed_size>();
    //!\brief Flag indicating whether debug mode is enabled.
    static constexpr bool is_debug = configuration_t::template exists<detail::debug_mode>();
    //!\brief Flag indicating whether a user provided callback was given.
    static constexpr bool is_one_way_execution = configuration_t::template exists<align_cfg::on_result>();
    //!\brief The selected scoring scheme.
    using scoring_scheme_type = decltype(get<align_cfg::scoring_scheme>(std::declval<configuration_t>()).value);
    //!\brief The alphabet of the selected scoring scheme.
    using scoring_scheme_alphabet_type = typename scoring_scheme_type::alphabet_type;
    //!\brief The selected result type.
    using result_type =
        std::remove_reference_t<decltype(seqan3::get<align_cfg::result>(std::declval<configuration_t>()))>;
    //!\brief The original score type selected by the user.
    using original_score_type = typename result_type::score_type;
    //!\brief The score type for the alignment algorithm.
    using score_type = std::conditional_t<is_vectorised, simd_type_t<original_score_type>, original_score_type>;
    //!\brief The trace directions type for the alignment algorithm.
    using trace_type = std::conditional_t<is_vectorised, simd_type_t<original_score_type>, trace_directions>;
    //!\brief The alignment result type if present. Otherwise seqan3::detail::empty_type.
    using alignment_result_type = decltype(determine_alignment_result_type());
    //!\brief The type of the matrix index.
    using matrix_index_type = std::conditional_t<is_vectorised,
                                                 simd_type_t<select_scalar_index_t<original_score_type>>,
                                                 size_t>;
    //!\brief The type of the matrix coordinate.
    using matrix_coordinate_type = lazy_conditional_t<is_vectorised,
                                                      lazy<simd_matrix_coordinate, matrix_index_type>,
                                                      matrix_coordinate>;

    //!\brief The number of alignments that can be computed in one simd vector.
    static constexpr size_t alignments_per_vector = [] () constexpr
                                                    {
                                                        if constexpr (is_vectorised)
                                                            return simd_traits<score_type>::length;
                                                        else
                                                            return 1;
                                                    }();
    //!\brief The rank of the selected result type.
    static constexpr int8_t result_type_rank = static_cast<int8_t>(decltype(std::declval<result_type>().value)::rank);
    //!\brief Flag indicating whether the score shall be computed.
    static constexpr bool compute_score = configuration_t::template exists<align_cfg::output_score_tag>();
    //!\brief Flag indicating whether the end positions shall be computed.
    static constexpr bool compute_end_positions =
        configuration_t::template exists<align_cfg::output_end_position_tag>();
    //!\brief Flag indicating whether the begin positions shall be computed.
    static constexpr bool compute_begin_positions =
        configuration_t::template exists<align_cfg::output_begin_position_tag>();
    //!\brief Flag indicating whether the sequence alignment shall be computed.
    static constexpr bool compute_sequence_alignment =
        configuration_t::template exists<align_cfg::output_alignment_tag>();
    //!\brief Flag indicating if any output option was set.
    static constexpr bool has_output_configuration = compute_score ||
                                                     compute_end_positions ||
                                                     compute_begin_positions ||
                                                     compute_sequence_alignment;
    //!\brief The padding symbol to use for the computation of the alignment.
    static constexpr original_score_type padding_symbol =
        static_cast<original_score_type>(1u << (sizeof_bits<original_score_type> - 1));
};

//------------------------------------------------------------------------------
// alignment_function_traits
//------------------------------------------------------------------------------

/*!\brief A traits class to provide a uniform access to the properties of the wrapped alignment algorithm.
 * \ingroup pairwise_alignment
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

}  // namespace seqan3::detail
