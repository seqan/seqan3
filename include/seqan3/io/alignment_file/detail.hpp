// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Auxiliary functions for the alignment IO.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <sstream>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/detail/to_string.hpp>
#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <range/v3/view/zip.hpp>

namespace seqan3::detail
{
//!\brief Comparator that is able two compare two views
struct view_equality_fn
{
    //!\brief Compares to ranges by delegating to std::ranges::equal.
    template <std::ranges::forward_range rng1_type, std::ranges::forward_range rng2_type>
    constexpr bool operator()(rng1_type && rng1, rng2_type && rng2) const
    {
        return std::ranges::equal(rng1, rng2);
    }
};

/*!\brief Compares two aligned sequence values and returns their cigar operation.
 * \ingroup alignment_file
 * \tparam reference_char_type Must be equality comparable to seqan3::gap.
 * \tparam query_char_type     Must be equality comparable to seqan3::gap.
 * \param  reference_char      The aligned character of the reference to compare.
 * \param  query_char          The aligned character of the query to compare.
 * \param  extended_cigar      Whether to print the extended cigar alphabet or not. See cigar operation.
 * \returns A seqan3::cigar_op representing the alignment operation between the two values.
 *
 * \details
 *
 * \note The resulting cigar operation is based on the query character (`query_char`).
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
 * In this case, the function seqan3::detail::map_aligned_values_to_cigar_op will return
 * 'D' since the query char is "deleted".
 *
 * The next alignment column shows the reference char ('C') on top and a
 * query char ('G') at the bottom.
 * ```
 * ... C ...
 *     |
 * ... G ...
 * ```
 * In this case, the function seqan3::detail::map_aligned_values_to_cigar_op will return
 * 'M', for the basic cigar the two bases are aligned, while
 * in the extended cigar alphabet (`extended_cigar` = `true`) the function
 * will return an 'X' since the bases are aligned but are not
 * equal.
 * \sa seqan3::aligned_sequence
 */
template <typename reference_char_type, typename query_char_type>
[[nodiscard]] constexpr cigar_op map_aligned_values_to_cigar_op(reference_char_type const reference_char,
                                                                query_char_type const query_char,
                                                                bool const extended_cigar)
//!\cond
    requires seqan3::detail::weakly_equality_comparable_with<reference_char_type, gap> &&
             seqan3::detail::weakly_equality_comparable_with<query_char_type, gap>
//!\endcond
{
    constexpr std::array<char, 6> operators{'M', 'D', 'I', 'P', 'X', '='};  // contains the possible cigar operators.
    uint8_t key = (static_cast<uint8_t>(reference_char == gap{}) << 1) | static_cast<uint8_t>(query_char == gap{});
    if (extended_cigar && (key == 0)) // in extended format refine the substitution operator to match/mismatch.
        key |= ((1 << 2) | static_cast<uint8_t>(query_char == reference_char));  // maps to [4, 5].

    return assign_char_to(operators[key], cigar_op{});
}

/*!\brief Creates a cigar string (SAM format) given an alignment represented by two aligned sequences.
 * \ingroup alignment_file
 *
 * \tparam alignment_type  Must model the seqan3::tuple_like and must have std::tuple_size 2.
 *                         Each tuple element must be a range over values comparable to seqan3::gap.
 * \param  alignment       The alignment, represented by a pair of aligned sequences,
 *                         to be transformed into cigar vector based on the
 *                         second (query) sequence.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See cigar operation.
 * \returns An std::vector<seqan3::cigar> representing the alignment.
 *
 * \details
 *
 * \note The resulting cigar_vector is based on the query sequence, which is the second sequence in the `alignment`
 *       pair.
 *
 * ### Example:
 *
 * Given the following alignment reference sequence on top and the query sequence at
 * the bottom:
 * ```
 * ATGG--CGTAGAGC
 * |||X  |||X|  |
 * ATGCCCCGTTG--C
 * ```
 * In this case, the function seqan3::detail::get_cigar_vector will return
 * the following cigar vector: "[('M',4),('I',2),('M',5),('D',2),('M',1)]".
 * The extended cigar string would look like this: "[('=',3)('X',1)('I',2)('=',3)('X',1)('=',1)('D',2)('=',1)]".
 *
 * \include test/snippet/io/alignment_file/get_cigar_vector.cpp
 *
 * \sa seqan3::aligned_sequence
 */
template <tuple_like alignment_type>
[[nodiscard]] inline std::vector<cigar> get_cigar_vector(alignment_type && alignment,
                                                         uint32_t const query_start_pos = 0,
                                                         uint32_t const query_end_pos = 0,
                                                         bool const extended_cigar = false)
//!\cond
    requires std::tuple_size_v<remove_cvref_t<alignment_type>> == 2
//!\endcond
{
    using std::get;

    auto & ref_seq = get<0>(alignment);
    auto & query_seq = get<1>(alignment);

    if (ref_seq.size() != query_seq.size())
        throw std::logic_error{"The aligned sequences must have the same length."};

    std::vector<cigar> result{};

    if (!ref_seq.size())
        return result; // return empty string if sequences are empty

    // Add (S)oft-clipping at the start of the read
    if (query_start_pos)
        result.emplace_back(query_start_pos, 'S'_cigar_op);

    // Create cigar string from alignment
    // -------------------------------------------------------------------------
    // initialize first operation and count value:
    cigar_op operation{map_aligned_values_to_cigar_op(ref_seq[0], query_seq[0], extended_cigar)};
    uint32_t count{0};

    // go through alignment columns
    for (auto column : views::zip(ref_seq, query_seq))
    {
        cigar_op next_op = map_aligned_values_to_cigar_op(std::get<0>(column), std::get<1>(column), extended_cigar);

        if (operation == next_op)
        {
            ++count;
        }
        else
        {
            result.emplace_back(count, operation);
            operation = next_op;
            count = 1;
        }
    }

    // append last cigar element
    result.emplace_back(count, operation);

    // Add (S)oft-clipping at the end of the read
    if (query_end_pos)
        result.emplace_back(query_end_pos, 'S'_cigar_op);

    return result;
}

/*!\brief Transforms a vector of cigar elements into a string representation.
 * \ingroup alignment_file
 * \param  cigar_vector The std::vector of seqan3::cigar elements to be transformed into a std::string.
 * \returns The cigar string (std::string).
 *
 * \details
 *
 * The transformation is done by printing the vector with the seqan3::debug_stream.
 */
[[nodiscard]] inline std::string get_cigar_string(std::vector<cigar> const & cigar_vector)
{
    std::string result{};
    std::ranges::for_each(cigar_vector, [&result] (auto & cig) { result.append(cig.to_string()); });
    return result;
}

/*!\brief Creates a cigar string (SAM format) given an alignment represented by two aligned sequences.
 * \ingroup alignment_file
 *
 * \tparam alignment_type  Must model the seqan3::tuple_like and must have std::tuple_size 2.
 *                         Each tuple element must be a range over values comparable to seqan3::gap.
 * \param  alignment       The alignment, represented by a pair of aligned sequences,
 *                         to be transformed into cigar vector based on the
 *                         second (query) sequence.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See cigar operation.
 * \returns An std::string representing the alignment as a cigar string.
 *
 * \details
 *
 * \note The resulting cigar_vector is based on the query sequence, which is the second sequence in the `alignment`
 *       pair.
 *
 * ### Example:
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
 * \sa seqan3::aligned_sequence
 */
template <tuple_like alignment_type>
[[nodiscard]] inline std::string get_cigar_string(alignment_type && alignment,
                                                  uint32_t const query_start_pos = 0,
                                                  uint32_t const query_end_pos = 0,
                                                  bool const extended_cigar = false)
//!\cond
    requires std::tuple_size_v<remove_cvref_t<alignment_type>> == 2
//!\endcond
{
    return get_cigar_string(get_cigar_vector(alignment, query_start_pos, query_end_pos, extended_cigar));
}

/*!\brief Transforms an alignment represented by two aligned sequences into the
 *        corresponding cigar string.
 * \ingroup alignment_file
 *
 * \tparam ref_seq_type    Must model std::ranges::forward_range. The value_type must
 *                         be equality comparable to seqan3::gap.
 * \tparam query_seq_type  Must model std::ranges::forward_range. The value_type must
 *                         be equality comparable to seqan3::gap.
 * \param  ref_seq         The reference sequence to compare against the query sequence.
 * \param  query_seq       The query sequence to build the cigar string for.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See cigar operation.
 * \returns An std::string representing the alignment as a cigar string.
 *
 * \details
 *
 * \note The resulting cigar string is based on the query sequence (`query_seq`).
 *
 * ### Example:
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
 *
 * \sa seqan3::aligned_sequence
 */
template <std::ranges::forward_range ref_seq_type, std::ranges::forward_range query_seq_type>
[[nodiscard]] inline std::string get_cigar_string(ref_seq_type && ref_seq,
                                                  query_seq_type && query_seq,
                                                  uint32_t const query_start_pos = 0,
                                                  uint32_t const query_end_pos = 0,
                                                  bool const extended_cigar = false)
//!\cond
    requires std::detail::weakly_equality_comparable_with<gap, std::ranges::range_reference_t<ref_seq_type>> &&
             std::detail::weakly_equality_comparable_with<gap, std::ranges::range_reference_t<query_seq_type>>
//!\endcond
{
    return get_cigar_string(std::tie(ref_seq, query_seq), query_start_pos, query_end_pos, extended_cigar);
}

/*!\brief Transforms a std::vector of operation-count pairs (representing the cigar string).
 * \ingroup alignment_file
 *
 * \tparam alignment_type The type of alignment; must model seqan3::tuple_like and all tuple element types
 *                        must model seqan3::aligned_sequence.
 *
 * \param[in,out] alignment    The alignment to fill with gaps according to the cigar information.
 * \param[in]     cigar_vector The cigar information given as a std::vector over seqan3::cigar.
 *
 * \details
 *
 * ### Example:
 *
 * Given the following cigar string "4M2I5M2D1M", the cigar information extracted by seqan3::detail::parse_cigar
 * would be "[(M,4), (I,2), (M,5), (D,2), (M,1)]". Given those cigar information, and an alignment variable containing
 * the two unaligned sequences "(ATGGCGTAGAGC, ATGCCCCGTTGC)", the alignment will be filled with the following gaps:
 *
 * ```
 * ATGG--CGTAGAGC
 * |||   ||| |  |
 * ATGCCCCGTTG--C
 * ```
 */
template <tuple_like alignment_type>
//!\cond
    requires std::tuple_size_v<remove_cvref_t<alignment_type>> == 2 &&
             detail::all_satisfy_aligned_seq<detail::tuple_type_list_t<alignment_type>>
//!\endcond
inline void alignment_from_cigar(alignment_type & alignment, std::vector<cigar> const & cigar_vector)
{
    using std::get;
    auto current_ref_pos  = std::ranges::begin(get<0>(alignment));
    auto current_read_pos = std::ranges::begin(get<1>(alignment));

    for (auto [cigar_count, cigar_operation] : cigar_vector)
    {
        // ignore since alignment shall contain sliced sequences
        if (('S'_cigar_op == cigar_operation) || ('H'_cigar_op == cigar_operation))
            continue;

        assert(('M'_cigar_op == cigar_operation) || ('='_cigar_op == cigar_operation) ||
               ('X'_cigar_op == cigar_operation) || ('D'_cigar_op == cigar_operation) ||
               ('N'_cigar_op == cigar_operation) || ('I'_cigar_op == cigar_operation) ||
               ('P'_cigar_op == cigar_operation)); // checked during IO

        if (('M'_cigar_op == cigar_operation) || ('='_cigar_op == cigar_operation) || ('X'_cigar_op == cigar_operation))
        {
            std::ranges::advance(current_ref_pos , cigar_count);
            std::ranges::advance(current_read_pos, cigar_count);
        }
        else if (('D'_cigar_op == cigar_operation) || ('N'_cigar_op == cigar_operation)) // insert gaps into read
        {
            assert(std::distance(current_read_pos, std::ranges::end(get<1>(alignment))) >= 0);
            current_read_pos = insert_gap(get<1>(alignment), current_read_pos, cigar_count);
            ++current_read_pos;
            std::ranges::advance(current_ref_pos , cigar_count);
        }
        else if (('I'_cigar_op == cigar_operation)) // Insert gaps into ref
        {
            assert(std::ranges::distance(current_ref_pos, std::ranges::end(get<0>(alignment))) >= 0);
            current_ref_pos = insert_gap(get<0>(alignment), current_ref_pos, cigar_count);
            ++current_ref_pos;
            std::ranges::advance(current_read_pos, cigar_count);
        }
        else if (('P'_cigar_op == cigar_operation)) // skip padding
        {
            current_ref_pos = insert_gap(get<0>(alignment), current_ref_pos, cigar_count);
            ++current_ref_pos;

            current_read_pos = insert_gap(get<1>(alignment), current_read_pos, cigar_count);
            ++current_read_pos;
        }
    }
}

//!\brief A functor that always throws when calling `operator()` (needed for the alignment "dummy" sequence).
struct access_restrictor_fn
{
    //!\brief Always throws a std::logic_error when called.
    template <typename chr_t>
    [[noreturn]] chr_t operator()(chr_t) const
    {
         throw std::logic_error{"Access is not allowed because there is no sequence information."};
    }
};

} // namespace seqan3::detail
