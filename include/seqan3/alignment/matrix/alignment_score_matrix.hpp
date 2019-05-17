// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
//!       declaration must satisfy seqan3::detail::Matrix.
//!\attention This is a pure base class, you must only use its specialisations.
//!\ingroup alignment_matrix
//!\implements seqan3::detail::Matrix
template <typename ...>
class alignment_score_matrix;

/*!\brief A score matrix represented in a one-dimensional std::vector
 * \ingroup alignment_matrix
 * \tparam score_t    The type of the score.
 *
 *
 * This data structure stores the matrix in a flat way using the
 * std::vector<#entry_type> data structure where each row is stored
 * continuously.
 *
 * # Example
 *
 * \snippet test/snippet/alignment/matrix/alignment_score_matrix.cpp code
 *
 * ### Output
 * \include test/snippet/alignment/matrix/alignment_score_matrix.out
 */
template <std::Integral score_t>
class alignment_score_matrix<std::vector<score_t>>
    : public row_wise_matrix<score_t>
{
public:
    using row_wise_matrix<score_t>::row_wise_matrix;
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_score_matrix
 * \{
 */

//!\brief Deduce the score matrix type from the provided arguments.
template <typename score_t>
alignment_score_matrix(std::vector<score_t>, size_t rows, size_t cols)
    -> alignment_score_matrix<std::vector<score_t>>;
//!\}

} // namespace seqan3::detail
