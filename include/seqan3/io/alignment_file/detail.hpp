
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
 * \brief Auxiliary functions for the alignment IO.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <sstream>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <range/v3/view/zip.hpp>

namespace seqan3::detail
{

/*!\brief Compares two aligned sequence values and returns their CIGAR operation.
 * \ingroup alignment_file
 * \relatesalso seqan3::aligned_sequence_concept
 * \tparam reference_char_type Must be equality comparable to seqan3::gap.
 * \tparam query_char_type     Must be equality comparable to seqan3::gap.
 * \param  reference_char      The aligned character of the reference to compare.
 * \param  query_char          The aligned character of the query to compare.
 * \param  extended_cigar      Whether to print the extended cigar alphabet or not. See CIGAR operation.
 * \returns A char representing the alignment operation between the
 *          two values.
 *
 * \details
 *
 * \attention Note that CIGAR elements (respectively by their
 *            CIGAR operation) are always related to on one of the two sequences
 *            in a pairwise alignment. In this case, the resulting
 *            operation is related to the \p query_char.
 *
 * ### Example:
 *
 * The following alignment column shows the reference char ('C') on top and a
 * gap for the query char at the bottom.
 * ```
 * ... C ...
 *     |
 * ... - ...
 * ```
 * In this case, the function seqan3::detail::compare_aligned_values will return
 * seqan3::'D', since the query char is "deleted".
 *
 * The next alignment column shows the reference char ('C') on top and a
 * query char ('G') at the bottom.
 * ```
 * ... C ...
 *     |
 * ... G ...
 * ```
 * In this case, the function seqan3::detail::compare_aligned_values will return
 * 'M', for the basic cigar the two bases are aligned, while
 * in the extended CIGAR alphabet (\p extended_cigar = `true`) the function
 * will return an 'X' since the bases are aligned but are not
 * equal.
 */
template<typename reference_char_type, typename query_char_type>
//!\cond
    requires std::EqualityComparableWith<std::remove_reference_t<reference_char_type>, gap> &&
             std::EqualityComparableWith<std::remove_reference_t<query_char_type>, gap>
//!\endcond
char compare_aligned_values(reference_char_type const reference_char,
                            query_char_type const query_char,
                            bool const extended_cigar)
{
    return (reference_char == gap::GAP)
                ? (query_char == gap::GAP)
                    ? 'P'
                    : 'I'
                : (query_char == gap::GAP)
                    ? 'D'
                    : (extended_cigar)
                        ? (query_char == reference_char)
                            ? '='
                            : 'X'
                        : 'M';
}

/*!\brief Transforms an alignment represented by two aligned sequences into the
 *        corresponding CIGAR string.
 * \ingroup alignment_file
 * \relatesalso seqan3::aligned_sequence_concept
 * \tparam ref_seq_type    Must model seqan3::aligned_sequence_concept.
 * \tparam query_seq_type  Must model seqan3::aligned_sequence_concept.
 * \param  alignment       The alignment, represented by a pair of aligned sequences,
 *                         to be transformed into CIGAR_vector based on the
 *                         second (query) sequence.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See CIGAR operation.
 * \returns An std::string representing the alignment.
 *
 * \details
 *
 * \attention Note that CIGAR elements (respectively by their
 *            CIGAR operation) are always related to on one of the two sequences
 *            in a pairwise alignment. In this case, the resulting cigar_vector
 *            is based on sequence at the second position of the \p alignment pair,
 *            namely the query sequence.
 *
 * ### Theoretical Example:
 *
 * The following alignment reference sequence on top and the query sequence at
 * the bottom.
 * ```
 * ATGG--CGTAGAGC
 * |||X  |||X|  |
 * ATGCCCCGTTG--C
 * ```
 * In this case, the function seqan3::detail::compare_aligned_values will return
 * The following cigar string when printed: "4M2I5M2D1M". The extended cigar
 * string would look like this: "3=1X2I3=1X1=2D1=".
 */
template<std::ranges::ForwardRange ref_seq_type, std::ranges::ForwardRange query_seq_type>
//\!cond
    requires std::EqualityComparableWith<gap, value_type_t<ref_seq_type>> &&
             std::EqualityComparableWith<gap, value_type_t<query_seq_type>>
//\!endcond
std::string get_cigar_string(std::pair<ref_seq_type, query_seq_type> const & alignment,
                             uint32_t query_start_pos = 0,
                             uint32_t query_end_pos = 0,
                             bool extended_cigar = false)
{
    if (std::get<0>(alignment).size() != std::get<1>(alignment).size())
        throw std::logic_error("The aligned sequences must have the same length.");

    std::ostringstream result;

    if (!std::get<0>(alignment).size())
        return std::string(); // return empty string if sequences are empty

    // Add (S)oft-clipping at the start of the read
    if (query_start_pos)
        result << query_start_pos << 'S';

    // Create cigar string from alignment
    // -------------------------------------------------------------------------
    // initialize first operation:
    char tmp_char{compare_aligned_values(std::get<0>(alignment)[0], std::get<1>(alignment)[0], extended_cigar)};
    unsigned tmp_length{0};

    // go through alignment columns
    for (auto column : ranges::view::zip(std::get<0>(alignment), std::get<1>(alignment)))
    {
        char next_op = compare_aligned_values(std::get<0>(column), std::get<1>(column), extended_cigar);

        if (tmp_char == next_op)
        {
            ++tmp_length;
        }
        else
        {
            result << tmp_length << tmp_char;
            tmp_char = next_op;
            tmp_length = 1;
        }
    }
    // append last cigar element
    result << tmp_length << tmp_char;

    // Add (S)oft-clipping at the end of the read
    if (query_end_pos)
        result << query_end_pos << 'S';

    return result.str();
}

} // namespace seqan3::detail
