// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Auxiliary functions for the alignment IO.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <sstream>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <range/v3/view/zip.hpp>

namespace seqan3::detail
{
/*!\brief Compares two aligned sequence values and returns their CIGAR operation.
 * \ingroup alignment_file
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
 * \sa seqan3::aligned_sequence_concept
 */
template<typename reference_char_type, typename query_char_type>
//!\cond
    requires std::detail::WeaklyEqualityComparableWith<reference_char_type, gap> &&
             std::detail::WeaklyEqualityComparableWith<query_char_type, gap>
//!\endcond
char compare_aligned_values(reference_char_type const reference_char,
                            query_char_type const query_char,
                            bool const extended_cigar)
{
    return (reference_char == gap{})
                ? (query_char == gap{})
                    ? 'P'
                    : 'I'
                : (query_char == gap{})
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
 *
 * \tparam ref_seq_type    Must model std::ranges::ForwardRange. The value_type must
 *                         be equality comparable to seqan3::gap.
 * \tparam query_seq_type  Must model std::ranges::ForwardRange. The value_type must
 *                         be equality comparable to seqan3::gap.
 * \param  ref_seq         The reference sequence to compare against the query sequence.
 * \param  query_seq       The query sequence to build the CIGAR string for.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See CIGAR operation.
 * \returns An std::string representing the alignment as a CIGAR string.
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
 * In this case, the function seqan3::detail::get_cigar_string will return
 * the following cigar string when printed: "4M2I5M2D1M". The extended cigar
 * string would look like this: "3=1X2I3=1X1=2D1=".
 * \sa seqan3::aligned_sequence_concept
 */
template<std::ranges::ForwardRange ref_seq_type, std::ranges::ForwardRange query_seq_type>
//!\cond
    requires std::detail::WeaklyEqualityComparableWith<gap, reference_t<ref_seq_type>> &&
             std::detail::WeaklyEqualityComparableWith<gap, reference_t<query_seq_type>>
//!\endcond
std::string get_cigar_string(ref_seq_type && ref_seq,
                             query_seq_type && query_seq,
                             uint32_t const query_start_pos = 0,
                             uint32_t const query_end_pos = 0,
                             bool const extended_cigar = false)
{
    if (ref_seq.size() != query_seq.size())
        throw std::logic_error("The aligned sequences must have the same length.");

    std::ostringstream result;

    if (!ref_seq.size())
        return std::string(); // return empty string if sequences are empty

    // Add (S)oft-clipping at the start of the read
    if (query_start_pos)
        result << query_start_pos << 'S';

    // Create cigar string from alignment
    // -------------------------------------------------------------------------
    // initialize first operation:
    char tmp_char{compare_aligned_values(ref_seq[0], query_seq[0], extended_cigar)};
    unsigned tmp_length{0};

    // go through alignment columns
    for (auto column : ranges::view::zip(ref_seq, query_seq))
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

/*!\brief Transforms an alignment represented by two aligned sequences into the
 *        corresponding CIGAR string.
 * \ingroup alignment_file
 *
 * \tparam alignment_type  Must model the seqan3::tuple_like_concept and must
 *                         have std::tuple_size 2. Each tuple element must model
 *                         std::ForwardRange and its value_type must be comparable
 *                         to seqan3::gap.
 * \param  alignment       The alignment, represented by a pair of aligned sequences,
 *                         to be transformed into CIGAR_vector based on the
 *                         second (query) sequence.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See CIGAR operation.
 * \returns An std::string representing the alignment as a CIGAR string.
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
 * In this case, the function seqan3::detail::get_cigar_string will return
 * the following cigar string when printed: "4M2I5M2D1M". The extended cigar
 * string would look like this: "3=1X2I3=1X1=2D1=".
 * \sa seqan3::aligned_sequence_concept
 */
template<tuple_like_concept alignment_type>
//!\cond
    requires std::tuple_size_v<remove_cvref_t<alignment_type>> == 2
//!\endcond
std::string get_cigar_string(alignment_type && alignment,
                             uint32_t const query_start_pos = 0,
                             uint32_t const query_end_pos = 0,
                             bool const extended_cigar = false)
{
    return get_cigar_string(get<0>(alignment), get<1>(alignment),
                            query_start_pos, query_end_pos, extended_cigar);
}

} // namespace seqan3::detail
