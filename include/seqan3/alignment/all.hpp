// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link alignment Alignment module \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup alignment Alignment
 * \brief The alignment module contains concepts, algorithms and classes that are related to the computation of
 *        pairwise and multiple sequence alignments.
 *
 * \details
 *
 * There are several types of alignments. We support \ref alignment_pairwise "pairwise alignments" so far, but we also
 * plan to support multiple alignments in the future.
 *
 * SeqAn offers a generic multi-purpose alignment library comprising all widely known alignment algorithms as well as
 * many special algorithms. These algorithms are all accessible through an easy to use alignment interface which
 * is described in \ref alignment_pairwise.
 */

#pragma once

#include <seqan3/alignment/aligned_sequence/all.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/decorator/all.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/matrix/all.hpp>
#include <seqan3/alignment/pairwise/all.hpp>
#include <seqan3/alignment/scoring/all.hpp>
