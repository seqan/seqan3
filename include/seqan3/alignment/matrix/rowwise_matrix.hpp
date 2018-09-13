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
 * \brief Contains seqan3::detail::rowwise_matrix.
 */

#pragma once

#include <vector>
#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief A matrix represented in a one-dimensional std::vector
 * \ingroup alignment_matrix
 * \implements seqan3::detail::matrix_concept
 * \tparam entry_t \copydoc seqan3::detail::rowwise_matrix::entry_type
 *
 * \details
 *
 * This data structure stores the matrix in a flat way using the
 * std::vector<#entry_type> data structure where each row is stored
 * continuously.
 */
template <typename entry_t>
struct rowwise_matrix
{
    //!\brief The type of the entry.
    using entry_type = entry_t;

    /*!\name Constructors, destructor and assignment
     * The copy-constructor, move-constructor, copy-assignment, move-assignment,
     * and destructor are implicitly defined.
     * \{
     */
     rowwise_matrix() = default;
     rowwise_matrix(rowwise_matrix const &) = default;
     rowwise_matrix(rowwise_matrix &&) = default;
     rowwise_matrix & operator=(rowwise_matrix const &) = default;
     rowwise_matrix & operator=(rowwise_matrix &&) = default;
    /*!\brief Construct the matrix out of the *entries*, the *rows*,
     *        and the *cols*.
     * \param entries The entry values as a flat std::vector <#entry_type>.
     * \param rows    The number of rows.
     * \param cols    The number of columns.
     */
    rowwise_matrix(std::vector<entry_type> entries, size_t rows, size_t cols)
        : _entries{std::move(entries)}, _rows{rows}, _cols{cols}
    {}
    //!\}

    //!\copydoc seqan3::detail::matrix_concept::rows
    inline std::size_t rows() const noexcept
    {
        return _rows;
    }

    //!\copydoc seqan3::detail::matrix_concept::cols
    inline std::size_t cols() const noexcept
    {
        return _cols;
    }

    //!\copydoc seqan3::detail::matrix_concept::at
    inline entry_type at(unsigned row, unsigned col) const noexcept
    {
        return _entries[row * cols() + col];
    }

private:
    //!\brief The matrix as a one-dimensional vector of entries
    //!       (each row is continuously stored).
    std::vector<entry_type> _entries;

    //!\copydoc seqan3::detail::matrix_concept::rows
    size_t _rows;

    //!\copydoc seqan3::detail::matrix_concept::cols
    size_t _cols;
};

} // namespace seqan3::detail
