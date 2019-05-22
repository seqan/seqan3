// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/view/take_until.hpp>
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
    template <std::ranges::ForwardRange rng1_type, std::ranges::ForwardRange rng2_type>
    constexpr bool operator()(rng1_type && rng1, rng2_type && rng2) const
    {
        return std::ranges::equal(rng1, rng2);
    }
};

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
 * 'D' since the query char is "deleted".
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
 * \sa seqan3::AlignedSequence
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
 * \sa seqan3::AlignedSequence
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
        throw std::logic_error{"The aligned sequences must have the same length."};

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
    for (auto column : std::view::zip(ref_seq, query_seq))
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

/*!\brief Creates a CIGAR string (SAM format) given an alignment represented by two aligned sequences.
 * \ingroup alignment_file
 *
 * \tparam alignment_type  Must model the seqan3::TupleLike and must
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
 * \sa seqan3::AlignedSequence
 */
template<TupleLike alignment_type>
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

/*!\brief Parses a cigar string into a vector of operation-count pairs (e.g. (M, 3)).
 * \ingroup alignment_file
 * \tparam cigar_input_type The type of a single pass input view over the cigar string; must model
 *                          std::ranges::InputRange.
 * \param[in]  cigar_input  The single pass input view over the cigar string to parse.
 *
 * \returns A tuple of size five containing (1) std::vector of operation-count pairs, e.g. (M, 3), that describe
 *          the alignment, (2) the aligned reference length, (3) the aligned query sequence length, (4) The number of
 *          soft clipped bases at the start of the alignment and (5) the number of bases at the end of the alignment.
 *
 * \details
 *
 * For example, the view over the cigar string "1S4M1D2M2S" will return `{[(M,4), (D,1), (M,2)], 7, 6, 1, 2}`.
 */
template <std::ranges::InputRange cigar_input_type>
std::tuple<std::vector<std::pair<char, size_t>>, size_t, size_t, size_t, size_t>
parse_cigar(cigar_input_type && cigar_input)
{
    std::vector<std::pair<char, size_t>> operations{};
    size_t sc_begin_count{}; // number of soft clipped bases at the beginning
    size_t sc_end_count{};   // number of soft clipped bases at the end
    std::array<char, 20> buffer{}; // buffer to parse numbers with from_chars. Biggest number should fit in uint64_t
    char cigar_op{'\0'};
    size_t cigar_count{0};
    size_t ref_length{0}, seq_length{0}; // length of aligned part for ref and query

    auto update_lengths_fn = [&ref_length, &seq_length, &cigar_op, &cigar_count] ()
        {
            if (is_char<'M'>(cigar_op) || is_char<'='>(cigar_op) || is_char<'X'>(cigar_op))
            {
                ref_length += cigar_count;
                seq_length += cigar_count;
            }
            else if (is_char<'D'>(cigar_op) || is_char<'N'>(cigar_op))
            {
                ref_length += cigar_count;
            }
            else if (is_char<'I'>(cigar_op))
            {
                seq_length += cigar_count;
            }
            else // illegal character
            {
                if (is_char<'P'>(cigar_op))
                    throw format_error{"We do currently not support cigar operation 'P'."};
                else
                    throw format_error{std::string{"Illegal cigar operation: "} + std::string{cigar_op}};
            }
        };

    // transform input into a signle input view if it isn't already
    auto cigar_view = cigar_input | view::single_pass_input;

    // check hard/soft clipping at the beginning manually
    // -----------------------------------------------------------------------------------------------------------------
    auto [ignore, buffer_end] = std::ranges::copy(cigar_view | view::take_until_or_throw(!is_digit), buffer.data());
    (void) ignore;

    cigar_op = *std::ranges::begin(cigar_view);
    std::ranges::next(std::ranges::begin(cigar_view));

    if (is_char<'H'>(cigar_op)) // hard clipping is ignored. parse the next operation
    {
        auto [ignore2, buffer_end2] = std::ranges::copy(cigar_view
                                                        | view::take_until_or_throw(!is_digit), buffer.data());
        buffer_end = buffer_end2;
        (void) ignore2;

        cigar_op = *std::ranges::begin(cigar_view);
        std::ranges::next(std::ranges::begin(cigar_view));
    }

    if (std::from_chars(buffer.begin(), buffer_end, cigar_count).ec != std::errc{})
        throw format_error{"Corrupted cigar string encountered"};

    if (is_char<'S'>(cigar_op)) // check for soft clipping at the beginning
    {
        sc_begin_count = cigar_count;
    }
    else
    {
        update_lengths_fn();
        operations.push_back({cigar_op, cigar_count});
    }

    // parse the rest of the cigar
    // -----------------------------------------------------------------------------------------------------------------
    while (std::ranges::begin(cigar_view) != std::ranges::end(cigar_view)) // until stream is not empty
    {
        buffer_end = (std::ranges::copy(cigar_view | view::take_until_or_throw(!is_digit), buffer.data())).out;
        cigar_op = *std::ranges::begin(cigar_view);
        std::ranges::next(std::ranges::begin(cigar_view));

        if (std::from_chars(buffer.begin(), buffer_end, cigar_count).ec != std::errc{})
            throw format_error{"Corrupted cigar string encountered"};

        if (is_char<'S'>(cigar_op)) // we are at the end, hard clipping afterwards can be ignored
        {
            sc_end_count = cigar_count;
            return {operations, ref_length, seq_length, sc_begin_count, sc_end_count};
        }
        update_lengths_fn();
        operations.push_back({cigar_op, cigar_count});
    }
    return {operations, ref_length, seq_length, sc_begin_count, sc_end_count};
}

/*!\brief Transforms a std::vector of operation-count pairs (representing the cigar string).
 * \ingroup alignment_file
 *
 * \tparam alignment_type The type of alignment; must model seqan3::TupleLike and all tuple element types
 *                        must model seqan3::AlignedSequence.
 *
 * \param[in,out] alignment  The alignment to fill with gaps according to the cigar information.
 * \param[in]     cigar      The cigar information given as a std::vector of operation-count pairs.
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
template <TupleLike alignment_type>
//!\cond
    requires std::tuple_size_v<remove_cvref_t<alignment_type>> == 2 &&
             detail::all_satisfy_aligned_seq<detail::tuple_type_list_t<alignment_type>>
//!\endcond
void alignment_from_cigar(alignment_type & alignment, std::vector<std::pair<char, size_t>> const & cigar)
{
    using std::get;
    auto current_ref_pos  = std::ranges::begin(get<0>(alignment));
    auto current_read_pos = std::ranges::begin(get<1>(alignment));

    for (auto [cigar_op, cigar_count] : cigar)
    {
        if (is_char<'M'>(cigar_op) || is_char<'='>(cigar_op) || is_char<'X'>(cigar_op))
        {
            std::ranges::advance(current_ref_pos , cigar_count);
            std::ranges::advance(current_read_pos, cigar_count);
        }
        else if (is_char<'D'>(cigar_op) || is_char<'N'>(cigar_op)) // insert gaps into read
        {
            assert(std::distance(current_read_pos, std::ranges::end(get<1>(alignment))) >= 0);
            current_read_pos = insert_gap(get<1>(alignment), current_read_pos, cigar_count);
            ++current_read_pos;
            std::ranges::advance(current_ref_pos , cigar_count);
        }
        else if (is_char<'I'>(cigar_op)) // Insert gaps into ref
        {
            assert(std::ranges::distance(current_ref_pos, std::ranges::end(get<0>(alignment))) >= 0);
            current_ref_pos = insert_gap(get<0>(alignment), current_ref_pos, cigar_count);
            ++current_ref_pos;
            std::ranges::advance(current_read_pos, cigar_count);
        }
        else // illegal character
        {
            if (is_char<'P'>(cigar_op))
                throw format_error{"We do currently not support cigar operation 'P'."};
            else
                throw format_error{std::string{"Illegal cigar operation: "} + std::string{cigar_op}};
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
