// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link alignment_matrix Matrix sub-module \endlink.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

/*!\defgroup alignment_matrix Matrix
 * \brief Provides data structures for representing alignment coordinates and alignments as a matrix.
 * \ingroup alignment
 *
 * \todo Write detailed landing page.
 */

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/row_wise_matrix.hpp>
