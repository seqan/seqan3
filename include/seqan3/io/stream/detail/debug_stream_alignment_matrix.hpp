// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides alignment matrices overloads for seqan3::debug_stream
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <sstream>

#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace seqan3::detail
{
template <Matrix alignment_matrix_t>
class alignment_matrix_formatter; // forward
} // namespace seqan3::detail

namespace seqan3
{

/*!\brief An alignment matrix can be printed to the seqan3::debug_stream.
 * \tparam alignment_matrix_t Type of the alignment matrix to be printed; must model seqan3::detail::Matrix.
 * \param s The seqan3::debug_stream.
 * \param matrix The alignment matrix.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * This prints out an alignment matrix, which can be a score matrix or a trace matrix.
 */
template <typename alignment_matrix_t>
//!\cond
    requires detail::Matrix<std::decay_t<alignment_matrix_t>>
//!\endcond
inline debug_stream_type & operator<<(debug_stream_type & s, alignment_matrix_t && matrix)
{
    std::string database(matrix.cols(), ' ');
    std::string query(matrix.rows(), ' ');

    std::stringstream sstream{};
    detail::alignment_matrix_formatter{matrix}.format(database, query, sstream, std::nullopt);
    s << sstream.str();
    return s;
}

} // namespace seqan3
