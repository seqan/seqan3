// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the function seqan3::cigar_from_alignment and a helper struct seqan3::cigar_clipped_bases.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3
{

/*!\brief Helper struct to specialise soft and hard clipping when using seqan3::cigar_from_alignment.
 * \ingroup cigar_conversion
 *
 * A CIGAR string might have hard or soft clipping at the front or back, e.g., `2H3S100M3S2H`.
 */
struct cigar_clipped_bases
{
    uint32_t hard_front{}; //!< The number of hard clipped bases at the front of the CIGAR string.
    uint32_t hard_back{};  //!< The number of hard clipped bases at the back of the CIGAR string.
    uint32_t soft_front{}; //!< The number of soft clipped bases at the front of the CIGAR string.
    uint32_t soft_back{};  //!< The number of soft clipped bases at the back of the CIGAR string.
};

/*!\brief Creates a CIGAR string (SAM format) given a seqan3::detail::pairwise_alignment represented by two
 *        `seqan3::aligned_sequence`s.
 * \ingroup cigar_conversion
 *
 * \tparam alignment_type  Must model seqan3::detail::pairwise_alignment.
 * \param  alignment       The alignment, represented by a pair of aligned sequences,
 *                         to be transformed into CIGAR vector based on the
 *                         second (query) sequence.
 * \param  clipped_bases Provides information on whether the query sequence was cropped (hard clipping) before the
 *                       alignment or whether part of the query sequence does not take part (soft clipping) in the
 *                       alignment.
 * \param  extended_cigar  Whether to print the extended CIGAR alphabet or not. See `seqan3::cigar::operation`.
 * \returns A std::vector\<seqan3::cigar\> representing the alignment.
 *
 * \details
 *
 * \note The resulting `cigar_vector` is based on the query sequence, which is the second sequence in the `alignment`
 *       pair.
 *
 * ### Example:
 *
 * Given the following alignment, reference sequence on top and the query sequence at
 * the bottom:
 * ```
 * ATGG--CGTAGAGCTT
 * |||X  |||X|  |||
 * ATGCCCCGTTG--CTT
 * ```
 * In this case, the function `seqan3::cigar_from_alignment` will return the following CIGAR vector:
 *
 * ```
 * [('M',4),('I',2),('M',5),('D',2),('M',3)]
 * ```
 *
 * The extended CIGAR string would look like this:
 *
 * ```
 * [('=',3)('X',1)('I',2)('=',3)('X',1)('=',1)('D',2)('=',3)]
 * ```
 *
 * \include test/snippet/alignment/cigar_conversion/cigar_from_alignment.cpp
 *
 * ### Soft and Hard clipping
 *
 * The terms soft and hard clipping were introduced by the
 * [SAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf). A SAM file is only storing a semi alignment
 * represented by the CIGAR string. The semi alignment of a query sequence is most often the result of a
 * [read mapping.](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/mapping/tutorial.html)
 *
 * #### Hard clipping
 *
 * Before aligning a query or read to the reference, some tools crop the query sequence because the quality is bad at
 * one end (e.g., Illumina reads tend to display a bad sequence quality towards the end of the read).
 *
 * To inform the user of a SAM file that query sequences were altered, hard-clipping information is appended to the
 * CIGAR string. E.g. `100M50H` indicates that a read of former length `150` has been cropped at the right end by `50`
 * bases. The sequence in the SAM file will thus only be 100 bases long.
 *
 * #### Soft clipping
 *
 * In contrast to hard clipping, soft clipping indicates that the read was cropped and the respective
 * bases do not participate in the alignment, but they are still part of the reported sequence.
 * E.g., `100M50S` indicates that a read of length `150` has been aligned without the rightmost `50` bases.
 * The sequence in the SAM file will still be 150 bases long.
 *
 * #### Adding soft and hard clipping with `seqan3::cigar_from_alignment`
 *
 * You can add the respective clipping information by passing an instance of seqan3::cigar_clipped_bases:
 *
 * \include test/snippet/alignment/cigar_conversion/cigar_from_alignment_with_clipping.cpp
 *
 * \sa seqan3::cigar_clipped_bases
 * \sa seqan3::sam_file_output
 * \sa seqan3::cigar
 * \sa seqan3::aligned_sequence
 */
template <typename alignment_type>
inline auto cigar_from_alignment(alignment_type const & alignment,
                                 cigar_clipped_bases const & clipped_bases = {},
                                 bool const extended_cigar = false)
{
    static_assert((tuple_like<std::remove_cvref_t<alignment_type>>
                   && std::tuple_size_v<std::remove_cvref_t<alignment_type>> == 2),
                  "The alignment must be a std::pair or std::tuple of size 2.");

    static_assert(
        (std::equality_comparable_with<gap, std::ranges::range_reference_t<decltype(std::get<0>(alignment))>>
         && std::equality_comparable_with<gap, std::ranges::range_reference_t<decltype(std::get<1>(alignment))>>),
        "The alignment must contain two ranges whose value_type is comparable to seqan3::gap.");

    /*brief: Compares two seqan3::aligned_sequence values and returns their cigar operation.
    * param  reference_char      The aligned character of the reference to compare.
    * param  query_char          The aligned character of the query to compare.
    * param  extended_cigar      Whether to print the extended cigar alphabet or not. See cigar operation.
    * returns A seqan3::cigar::operation representing the alignment operation between the two values.
    *
    * The resulting CIGAR operation is based on the query character (`query_char`).
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
    * In this case, `to_cigar_op` will return
    * 'D' since the query char is "deleted".
    *
    * The next alignment column shows the reference char ('C') on top and a
    * query char ('G') at the bottom.
    * ```
    * ... C ...
    *     |
    * ... G ...
    * ```
    * In this case, `to_cigar_op` will return
    * 'M', for the basic CIGAR the two bases are aligned, while
    * in the extended CIGAR alphabet (`extended_cigar` = `true`) the function
    * will return an 'X' since the bases are aligned but are not
    * equal.
    * \sa seqan3::aligned_sequence
    */
    auto to_cigar_op = [extended_cigar](auto const reference_char, auto const query_char)
    {
        // note that N is not considered because it is equivalent to D but has a special meaning:
        // SAM spec: "For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments,
        //            the interpretation of N is not defined."
        // as we cannot know the meaning, the user has to change D -> N themself
        constexpr std::array<char, 6> operators{'M', 'D', 'I', 'P', 'X', '='}; // contains the possible cigar operators.

        // no gaps               -> 00 -> 0
        // gap only in query     -> 01 -> 1
        // gap only in reference -> 10 -> 2
        // both gaps             -> 11 -> 3
        uint8_t key = (static_cast<uint8_t>(reference_char == gap{}) << 1) | static_cast<uint8_t>(query_char == gap{});

        // mismatch -> 100 -> 4
        // match    -> 101 -> 5
        if (extended_cigar && (key == 0)) // in extended format, refine the substitution operator to match/mismatch.
            key |= ((1 << 2) | static_cast<uint8_t>(query_char == reference_char)); // maps to [4, 5].

        return assign_char_to(operators[key], cigar::operation{});
    };

    using std::get;

    auto & ref_seq = get<0>(alignment);
    auto & query_seq = get<1>(alignment);

    if (ref_seq.size() != query_seq.size())
        throw std::logic_error{"The aligned sequences (including gaps) must have the same length."};

    if (std::ranges::empty(ref_seq)) // only check ref_seq because query_seq was checked to have to same size
        throw std::logic_error{"The aligned sequences may not be empty."};

    std::vector<cigar> result{};

    // Add (H)ard-clipping at the start of the query
    if (clipped_bases.hard_front)
        result.emplace_back(clipped_bases.hard_front, 'H'_cigar_operation);

    // Add (S)oft-clipping at the start of the query
    if (clipped_bases.soft_front)
        result.emplace_back(clipped_bases.soft_front, 'S'_cigar_operation);

    // Create cigar string from alignment
    // -------------------------------------------------------------------------
    // initialize first operation and count value:
    cigar::operation operation{to_cigar_op(ref_seq[0], query_seq[0])};
    uint32_t count{0};

    // go through alignment columns
    for (auto && [reference_char, query_char] : views::zip(ref_seq, query_seq))
    {
        cigar::operation next_op = to_cigar_op(reference_char, query_char);

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

    // Add (S)oft-clipping at the end of the query
    if (clipped_bases.soft_back)
        result.emplace_back(clipped_bases.soft_back, 'S'_cigar_operation);

    // Add (H)ard-clipping at the end of the query
    if (clipped_bases.hard_back)
        result.emplace_back(clipped_bases.hard_back, 'H'_cigar_operation);

    return result;
}

} // namespace seqan3
