// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link alignment_matrix matrix sub-module \endlink.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

/*!\defgroup alignment_matrix Matrix
 * \brief Provides data structures for representing alignment coordinates and alignments as a matrix.
 * \ingroup alignment
 * \see alignment
 */

#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/alignment/matrix/debug_matrix.hpp>
#ifdef SEQAN3_DEPRECATED_310
// for deprecated seqan3::alignment_coordinate
#include <seqan3/alignment/matrix/detail/advanceable_alignment_coordinate.hpp>
#endif // SEQAN3_DEPRECATED_310
#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/row_wise_matrix.hpp>
