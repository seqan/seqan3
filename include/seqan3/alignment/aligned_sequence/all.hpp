// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link alignment_aligned_sequence Alignment / Aligned Sequence submodule \endlink.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

/*!\defgroup alignment_aligned_sequence Aligned Sequence
 * \ingroup alignment
 * \brief Provides seqan3::aligned_sequence, as well as various ranges that model it.
 * \see alignment
 *
 * The seqan3::aligned_sequence concept can be used to describe a sequence that is augmented with gaps, e.g.
 * a sequence that is part of an alignment.
 *
 * The data structure that we use most often to model `seqan3::aligned_sequence` is the `seqan3::gap_decorator`.
 * It is a lightweight data structure that only holds a view on the sequence (no copy is made) and on top can hold
 * `seqan3::gap`s.
 *
 * E.g. `AC-GA` is a aligned sequences where `-` represents a gap.
 *
 * \sa alignment_pairwise
 */

#pragma once

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/aligned_sequence/debug_stream_alignment.hpp>
