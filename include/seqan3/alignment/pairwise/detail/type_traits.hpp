// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides helper type traits for the configuration and execution of the alignment algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

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

/*!\brief A traits type for the alignment algorithm that exposes static information stored within the alignment
 *        configuration object.
 * \ingroup pairwise_alignment
 *
 * \tparam config_t The type of the alignment configuration object; must be a specialisation of seqan3::configuration.
 */
template <typename config_t>
//!\cond
    requires is_type_specialisation_of_v<config_t, configuration>
//!\endcond
struct alignment_configuration_traits
{
    //!\brief Flag to indicate vectorised mode.
    static constexpr bool is_vectorised = config_t::template exists<remove_cvref_t<decltype(align_cfg::vectorise)>>();
    //!\brief Flag indicating whether parallel alignment mode is enabled.
    static constexpr bool is_parallel = config_t::template exists<align_cfg::parallel>();
    //!\brief Flag indicating whether local alignment mode is enabled.
    static constexpr bool is_local = config_t::template exists<align_cfg::mode<detail::local_alignment_type>>();
    //!\brief Flag indicating whether banded alignment mode is enabled.
    static constexpr bool is_banded = config_t::template exists<align_cfg::band>();
    //!\brief Flag indicating whether debug mode is enabled.
    static constexpr bool is_debug = config_t::template exists<detail::debug_mode>();

    //!\brief The configured alignment mode.
    using alignment_mode_t = decltype(get<align_cfg::mode>(std::declval<config_t>()).value);
    //!\brief The selected scoring scheme.
    using scoring_scheme_t = decltype(get<align_cfg::scoring>(std::declval<config_t>()).value);
    //!\brief The alphabet of the selected scoring scheme.
    using scoring_scheme_alphabet_t = typename scoring_scheme_t::alphabet_type;
    //!\brief The selected result type.
    using result_t = std::remove_reference_t<decltype(seqan3::get<align_cfg::result>(std::declval<config_t>()))>;
    //!\brief The original score type selected by the user.
    using original_score_t = typename result_t::score_type;
    //!\brief The score type for the alignment algorithm.
    using score_t = std::conditional_t<is_vectorised, simd_type_t<original_score_t>, original_score_t>;
    //!\brief The trace directions type for the alignment algorithm.
    using trace_t = std::conditional_t<is_vectorised, simd_type_t<original_score_t>, trace_directions>;

    //!\brief The number of alignments that can be computed in one simd vector.
    static constexpr size_t alignments_per_vector = [] () constexpr
                                                    {
                                                        if constexpr (is_vectorised)
                                                            return simd_traits<score_t>::length;
                                                        else
                                                            return 1;
                                                    }();
    //!\brief The rank of the selected result type.
    static constexpr int8_t result_type_rank = static_cast<int8_t>(decltype(std::declval<result_t>().value)::rank);
};

}  // namespace seqan3::detail
