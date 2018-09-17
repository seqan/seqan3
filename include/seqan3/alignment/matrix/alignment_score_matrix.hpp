// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains the declaration of seqan3::detail::alignment_score_matrix.
 */

#pragma once

#include <seqan3/alignment/matrix/row_wise_matrix.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

//!\brief The declaration of alignment_score_matrix. Each specialisation of this
//!       declaration must satisfy seqan3::detail::matrix_concept.
//!\attention This is a pure base class, you must only use its specialisations.
//!\ingroup alignment_matrix
//!\implements seqan3::detail::matrix_concept
template <typename ...>
struct alignment_score_matrix;

/*!\brief A score matrix represented in a one-dimensional std::vector
 * \ingroup alignment_matrix
 * \tparam score_t    The type of the score.
 *
 *
 * This data structure stores the matrix in a flat way using the
 * std::vector<#entry_type> data structure where each row is stored
 * continuously.
 *
 * ## Example
 *
 * \snippet test/snippet/alignment/matrix/alignment_score_matrix.cpp code
 *
 * ### Output
 * \include test/snippet/alignment/matrix/alignment_score_matrix.out
 */
template <std::Integral score_t>
struct alignment_score_matrix<std::vector<score_t>>
    : public row_wise_matrix<score_t>
{
    using row_wise_matrix<score_t>::row_wise_matrix;
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_score_matrix
 * \{
 */
template <typename score_t>
alignment_score_matrix(std::vector<score_t>, size_t rows, size_t cols)
    -> alignment_score_matrix<std::vector<score_t>>;
//!\}

} // namespace seqan3::detail
