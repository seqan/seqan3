// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::row_wise_matrix.
 */

#pragma once

#include <vector>
#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief A matrix represented in a one-dimensional std::vector
 * \ingroup alignment_matrix
 * \implements seqan3::detail::Matrix
 * \tparam entry_t \copydoc seqan3::detail::row_wise_matrix::entry_type
 *
 * \details
 *
 * This data structure stores the matrix in a flat way using the
 * std::vector<#entry_type> data structure where each row is stored
 * continuously.
 */
template <typename entry_t>
class row_wise_matrix
{
public:
    //!\brief The type of the entry.
    using entry_type = entry_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    row_wise_matrix() = default;                                    //!< Defaulted
    row_wise_matrix(row_wise_matrix const &) = default;             //!< Defaulted
    row_wise_matrix(row_wise_matrix &&) = default;                  //!< Defaulted
    row_wise_matrix & operator=(row_wise_matrix const &) = default; //!< Defaulted
    row_wise_matrix & operator=(row_wise_matrix &&) = default;      //!< Defaulted
    ~row_wise_matrix() = default;                                   //!< Defaulted
    /*!\brief Construct the matrix out of the *entries*, the *rows*,
     *        and the *cols*.
     * \param entries The entry values as a flat std::vector <#entry_type>.
     * \param rows    The number of rows.
     * \param cols    The number of columns.
     */
    row_wise_matrix(std::vector<entry_type> entries, size_t const rows, size_t const cols)
        : _entries{std::move(entries)}, _rows{rows}, _cols{cols}
    {}
    //!\}

    //!\copydoc seqan3::detail::Matrix::rows
    size_t rows() const noexcept
    {
        return _rows;
    }

    //!\copydoc seqan3::detail::Matrix::cols
    size_t cols() const noexcept
    {
        return _cols;
    }

    //!\copydoc seqan3::detail::Matrix::at
    entry_type at(size_t const row, size_t const col) const noexcept
    {
        assert(row < rows() && col < cols());
        return _entries[row * cols() + col];
    }

private:
    //!\brief The matrix as a one-dimensional vector of entries
    //!       (each row is continuously stored).
    std::vector<entry_type> _entries;

    //!\copydoc seqan3::detail::Matrix::rows
    size_t _rows;

    //!\copydoc seqan3::detail::Matrix::cols
    size_t _cols;
};

} // namespace seqan3::detail
