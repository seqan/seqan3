// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Auxiliary functions for the SAM IO.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <ranges>
#include <seqan3/std/charconv>
#include <sstream>

#include <seqan3/alignment/detail/pairwise_alignment_concept.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{
//!\brief Comparator that is able two compare two views
//!\ingroup io_sam_file
struct view_equality_fn
{
    //!\brief Compares to ranges by delegating to std::ranges::equal.
    template <std::ranges::forward_range rng1_type, std::ranges::forward_range rng2_type>
    constexpr bool operator()(rng1_type && rng1, rng2_type && rng2) const
    {
        return std::ranges::equal(rng1, rng2);
    }
};

/*!\brief Updates the sequence lengths by `cigar_count` depending on the cigar operation `op`.
 * \ingroup io_sam_file
 * \param[in, out]  ref_length      The reference sequence's length.
 * \param[in, out]  seq_length      The query sequence's length.
 * \param[in]       cigar_operation The cigar operation.
 * \param[in]       cigar_count     The cigar count value to add to the length depending on the cigar operation.
 */
inline void update_alignment_lengths(int32_t & ref_length,
                                     int32_t & seq_length,
                                     char const cigar_operation,
                                     uint32_t const cigar_count)
{
    switch (cigar_operation)
    {
    case 'M':
    case '=':
    case 'X':
        ref_length += cigar_count, seq_length += cigar_count;
        break;
    case 'D':
    case 'N':
        ref_length += cigar_count;
        break;
    case 'I':
        seq_length += cigar_count;
        break;
    case 'S':
    case 'H':
    case 'P':
        break; // no op (soft-clipping or padding does not increase either length)
    default:
        throw format_error{"Illegal cigar operation: " + std::string{cigar_operation}};
    }
}

/*!\brief Parses a cigar string into a vector of operation-count pairs (e.g. (M, 3)).
 * \ingroup io_sam_file
 * \param[in]  cigar_str  The cigar string to parse.
 *
 * \returns A std::vector over seqan3::cigar that describes the alignment.
 *
 * \details
 *
 * For example, the view over the cigar string "1H4M1D2M2S" will return
 * `{[(H,1), (M,4), (D,1), (M,2), (S,2)], 7, 6}`.
 */
SEQAN3_WORKAROUND_LITERAL std::vector<cigar> parse_cigar(std::string_view const cigar_str)
{
    std::vector<seqan3::cigar> cigar_vector{};

    if (cigar_str == "*")
        return cigar_vector;

    uint32_t cigar_count{};
    char const * ptr = cigar_str.data();
    char const * const end = ptr + cigar_str.size();

    while (ptr < end)
    {
        auto const res = std::from_chars(ptr, end, cigar_count); // reads number up to next character

        if (res.ec != std::errc{})
            throw format_error{"Corrupted cigar string."};

        ptr = res.ptr + 1; // skip cigar operation character

        cigar_vector.emplace_back(cigar_count, seqan3::assign_char_strictly_to(*res.ptr, seqan3::cigar::operation{}));
    }

    return cigar_vector;
}

/*!\brief Transforms a vector of cigar elements into a string representation.
 * \ingroup io_sam_file
 * \param  cigar_vector The std::vector of seqan3::cigar elements to be transformed into a std::string.
 * \returns The cigar string (std::string).
 */
[[nodiscard]] inline std::string get_cigar_string(std::vector<cigar> const & cigar_vector)
{
    std::string result{};
    std::ranges::for_each(cigar_vector,
                          [&result](auto & cig)
                          {
                              result.append(static_cast<std::string_view>(cig.to_string()));
                          });
    return result;
}

/*!\brief Transforms an alignment represented by two seqan3::aligned_sequence's into the corresponding cigar string.
 * \ingroup io_sam_file
 *
 * \tparam ref_seq_type    Must model seqan3::aligned_sequence.
 * \tparam query_seq_type  Must model seqan3::aligned_sequence.
 * \param  ref_seq         The reference sequence to compare against the query sequence.
 * \param  query_seq       The query sequence to build the cigar string for.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See cigar operation.
 * \returns A std::string representing the alignment as a cigar string.
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
template <seqan3::aligned_sequence ref_seq_type, seqan3::aligned_sequence query_seq_type>
[[nodiscard]] inline std::string get_cigar_string(ref_seq_type && ref_seq,
                                                  query_seq_type && query_seq,
                                                  uint32_t const query_start_pos = 0,
                                                  uint32_t const query_end_pos = 0,
                                                  bool const extended_cigar = false)
{
    return get_cigar_string(std::tie(ref_seq, query_seq), query_start_pos, query_end_pos, extended_cigar);
}

//!\brief A functor that always throws when calling `operator()` (needed for the alignment "dummy" sequence).
//!\ingroup io_sam_file
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
