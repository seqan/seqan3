// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configuration for alignment output.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/concept/core_language.hpp>

namespace seqan3::align_cfg
{

/*!\brief Tag representing score output for the alignment algorithms.
 * \ingroup alignment_configuration
 */
struct output_score_tag : public pipeable_config_element<output_score_tag>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_score};
};

/*!\brief Configures the alignment result to output the score.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to compute and output the score. If this option is not set in the alignment
 * configuration, accessing the score via the seqan3::alignment_result object is forbidden and will lead to a compile
 * time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_score.cpp
 *
 * \see seqan3::align_cfg::output_end_position
 * \see seqan3::align_cfg::output_begin_position
 * \see seqan3::align_cfg::output_alignment
 * \see seqan3::align_cfg::output_sequence1_id
 * \see seqan3::align_cfg::output_sequence2_id
 */
inline constexpr output_score_tag output_score{};

/*!\brief Tag representing end position output for the alignment algorithms.
 * \ingroup alignment_configuration
 */
struct output_end_position_tag : public pipeable_config_element<output_end_position_tag>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_end_position};
};

/*!\brief Configures the alignment result to output the end position.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to compute and output the end positions of the aligned sequences.
 * The end positions must not be identical to the end of the original source sequences. For example,
 * the optimal local alignment might only represent a slice of the original sequences.
 * The end positions denote the end of the alignment within the original sequences, i.e. the positions behind
 * the last aligned characters.
 *
 * If this option is not set in the alignment configuration, then accessing the end positions via the
 * seqan3::alignment_result object is forbidden and will lead to a compile time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_end_position.cpp
 *
 * \see seqan3::align_cfg::output_score
 * \see seqan3::align_cfg::output_begin_position
 * \see seqan3::align_cfg::output_alignment
 * \see seqan3::align_cfg::output_sequence1_id
 * \see seqan3::align_cfg::output_sequence2_id
 */
inline constexpr output_end_position_tag output_end_position{};

} // namespace seqan3::align_cfg
