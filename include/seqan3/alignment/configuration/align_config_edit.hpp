// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_cfg::edit_scheme.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/hamming_scoring_scheme.hpp>
#include <seqan3/core/configuration/configuration.hpp>

namespace seqan3::align_cfg
{

/*!\brief Shortcut for edit distance configuration.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * The edit distance computation is a specific sub-problem of the alignment computation with the aim to count the number
 * of edits to transform one sequence into another. An edit operation can be a substitution, an insertion, or a
 * deletion. Accordingly, this algorithm uses a predefined scoring scheme as well as a gap scheme, where the score for
 * a match is 0, for a mismatch -1, for a gap -1, and for a gap open 0.
 *
 * ### Performance
 *
 * Under the hood SeqAn uses a [fast bit-vector algorithm](https://doi.org/10.1145/316542.316550) to compute the
 * edit distance whenever possible. This depends on the final alignment configuration. Currently, the fast
 * edit distance algorithm is only triggered for \ref seqan3::align_cfg::method_global "global alignments" with the
 * with free ends in the first sequence. So make sure to configure the seqan3::align_cfg::method_global configuration
 * element accordingly (see class documentation).
 *
 * The performance of the algorithm can further be improved if the number of maximal errors (edits) is known by using
 * the align_cfg::min_score configuration.
 *
 * \include snippet/alignment/configuration/align_cfg_edit_example.cpp
 *
 * \attention If the edit distance configuration is combined with any other configuration element or setting, the
 * algorithm falls back to the slower standard pairwise algorithm. For example the `cfg_slow` in the above example will
 * trigger the slower algorithm which can handle the case if the ends are free in the second sequence instead of the
 * first sequence.
 */
inline constexpr configuration edit_scheme =
    scoring_scheme{hamming_scoring_scheme{}} | gap_cost_affine{open_score{0}, extension_score{-1}};

} // namespace seqan3::align_cfg
