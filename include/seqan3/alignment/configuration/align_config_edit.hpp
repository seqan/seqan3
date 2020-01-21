// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::edit.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

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
 * Under the hood SeqAn uses a [fast bit-vector algorithm](https://dl.acm.org/citation.cfm?id=316550) to compute the
 * edit distance whenever possible. This depends on the final alignment configuration. Currently, the fast
 * edit distance algorithm is only triggered for \ref seqan3::global_alignment "global alignments" and
 * \ref seqan3::end_gaps::free_ends_first "semi-global alignments" with free ends in the first sequence.
 *
 * The performance of the algorithm can further be improved if the number of maximal errors (edits) is known by using
 * the align_cfg::max_error configuration.
 *
 * \include snippet/alignment/configuration/align_cfg_edit_example.cpp
 *
 * \attention If the edit distance configuration is combined with any other configuration element or setting, the
 * algorithm falls back to the slower standard pairwise algorithm. For example the `cfg_slow` in the above example will
 * trigger the slower algorithm which can handle the case if the ends are free in the second sequence instead of the
 * first sequence.
 */
inline constexpr configuration edit = mode{global_alignment} | scoring{nucleotide_scoring_scheme{}} |
                                      gap{gap_scheme{gap_score{-1}}};

} // namespace seqan3
