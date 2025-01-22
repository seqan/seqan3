// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the function seqan3::alignment_from_cigar.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/decorator/gap_decorator.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/utility/views/slice.hpp>

namespace seqan3
{

/*!\brief Construct an alignment from a CIGAR string and the corresponding sequences.
 * \ingroup cigar_conversion
 * \tparam reference_type The type of the reference sequence for a SAM record.
 * \tparam sequence_type The type of the read sequence for a SAM record.
 * \param[in] cigar_vector The CIGAR information to convert to an alignment.
 * \param[in] reference The reference sequence to which the `query` was aligned to, the alignment being represented by `cigar_vector`.
 * \param[in] zero_based_reference_start_position The start position of the alignment in the reference sequence. The
 *                                                position is zero-based. When using our seqan3::sam_file_input, note
 *                                                that the seqan3::sam_file_input::record_type::reference_position()
 *                                                is always zero-based.
 * \param[in] query The query or read sequence of the alignment represented by `cigar_vector`.
 * \returns An alignment represented by a `std::tuple` of size 2 holding 2 `seqan3::gap_decorator`s. At position 0 is
 *          the aligned reference sequence and at position 1 the aligned read sequence.
 *
 * ### Quick background on the CIGAR string
 *
 * The CIGAR string is a compact representation of an aligned read against a reference and was introduced by
 * the [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format. The SAM format stores the result of mapping
 * short/long read sequences from a sequencing experiment (e.g., Illumina/Nanopore) against a reference (e.g., hg38).
 *
 * ### Conversion to an alignment
 *
 * You can reconstruct a full alignment from a CIGAR string, if you have the respective sequences at hand:
 *
 * \include test/snippet/alignment/cigar_conversion/alignment_from_cigar.cpp
 *
 * ### Quick explanation of the alignment representation
 *
 * In seqan3, an alignment is represented by a `std::tuple` of size 2 that holds 2 `seqan3::aligned_sequence`s.
 *
 * The data structure that we use most often to model `seqan3::aligned_sequence` is the `seqan3::gap_decorator`.
 * It is a lightweight data structure that only holds a view on the sequence (no copy is made) and on top can hold
 * `seqan3::gap`s.
 *
 * In the above example, the read sequence `ACGT` is aligned to the reference with one gap:
 * `AC-GA` where `-` represents a gap.
 * In the CIGAR string, the gap in the query/read is represented by `1D`.
 *
 * The full alignment consist of two aligned sequences (read and reference).
 * In the above example, the alignment
 * ```
 * position   01234
 * reference  ACTGA
 * read       AC-GA
 * ```
 * is represented by a tuple of the aligned reference at the first position (std::get<0>) and the aligned read at the
 * second position (std::get<1>): `(ACTGA,AC-GA)`.
 *
 * ### IO Example
 *
 * A more realistic example is extracting the information directly from a SAM file:
 *
 * \include test/snippet/alignment/cigar_conversion/alignment_from_cigar_io.cpp
 *
 * \sa seqan3::sam_file_input
 * \sa seqan3::cigar
 */
template <typename reference_type, typename sequence_type>
inline auto alignment_from_cigar(std::vector<cigar> const & cigar_vector,
                                 reference_type const & reference,
                                 uint32_t const zero_based_reference_start_position,
                                 sequence_type const & query)
{
    if (cigar_vector.empty())
        throw std::logic_error{"An empty CIGAR is not a valid alignment representation."};

    // compute the length of the aligned region in the reference sequence
    // -------------------------------------------------------------------------
    // this requires a first stream over the cigar vector.
    uint32_t reference_length{0};
    uint32_t query_length{0};

    for (auto [cigar_count, cigar_operation] : cigar_vector)
    {
        if (('M'_cigar_operation == cigar_operation) || ('='_cigar_operation == cigar_operation)
            || ('X'_cigar_operation == cigar_operation))
        {
            reference_length += cigar_count;
            query_length += cigar_count;
        }
        else if ('D'_cigar_operation == cigar_operation)
        {
            reference_length += cigar_count;
        }
        else if ('I'_cigar_operation == cigar_operation)
        {
            query_length += cigar_count;
        }
    }

    if (static_cast<size_t>(zero_based_reference_start_position + reference_length) > std::ranges::size(reference))
        throw std::logic_error{"The CIGAR string indicates a reference length of at least "
                               + std::to_string(zero_based_reference_start_position + reference_length)
                               + ", but the supplied reference sequence is only of size"
                               + std::to_string(std::ranges::size(reference)) + "."};

    // Get soft clipping at the start and end of the CIGAR string
    // -------------------------------------------------------------------------
    uint32_t soft_clipping_start{0};
    uint32_t soft_clipping_end{0};

    // Checks whether the given index in the cigar vector is a soft clip.
    auto soft_clipping_at = [&cigar_vector](size_t const index)
    {
        return cigar_vector[index] == 'S'_cigar_operation;
    };
    // Checks whether the given index in the cigar vector is a hard clip.
    auto hard_clipping_at = [&](size_t const index)
    {
        return cigar_vector[index] == 'H'_cigar_operation;
    };
    // Checks whether the given cigar vector has at least min_size many elements.
    auto vector_size_at_least = [&](size_t const min_size)
    {
        return cigar_vector.size() >= min_size;
    };
    // Returns the cigar count of the ith cigar element in the given cigar vector.
    auto cigar_count_at = [&](size_t const index)
    {
        return get<0>(cigar_vector[index]);
    };

    // check for soft clipping at the first two positions
    // cigar is non-empty, checked at the very beginning.
    if (soft_clipping_at(0))
        soft_clipping_start = cigar_count_at(0);
    else if (vector_size_at_least(2) && hard_clipping_at(0) && soft_clipping_at(1))
        soft_clipping_start = cigar_count_at(1);

    // Check for soft clipping at the last two positions to validate the CIGAR string.
    // Even if the two following arithmetics overflow, they are protected by the corresponding if expressions below.
    auto const last_index = cigar_vector.size() - 1;
    auto const second_last_index = last_index - 1;

    if (vector_size_at_least(2) && soft_clipping_at(last_index))
        soft_clipping_end = cigar_count_at(last_index);
    else if (vector_size_at_least(3) && hard_clipping_at(last_index) && soft_clipping_at(second_last_index))
        soft_clipping_end = cigar_count_at(second_last_index);

    if (soft_clipping_start + query_length + soft_clipping_end != std::ranges::size(query))
        throw std::logic_error{"The CIGAR string indicates a query/read sequence length of "
                               + std::to_string(soft_clipping_start + query_length + soft_clipping_end)
                               + ", but the supplied query/read sequence is of size"
                               + std::to_string(std::ranges::size(query)) + "."};

    // Assign the sequence to the alignment (a tuple of 2 gap decorators)
    // -------------------------------------------------------------------------
    using gapped_reference_type = gap_decorator<decltype(reference | views::slice(0, 0))>;
    using gapped_sequence_type = gap_decorator<decltype(query | views::slice(0, 0))>;
    using alignment_type = std::tuple<gapped_reference_type, gapped_sequence_type>;

    alignment_type alignment{};

    assign_unaligned(get<0>(alignment),
                     reference
                         | views::slice(zero_based_reference_start_position,
                                        zero_based_reference_start_position + reference_length));
    // query_length already accounts for soft clipping at begin and end
    assign_unaligned(get<1>(alignment), query | views::slice(soft_clipping_start, soft_clipping_start + query_length));

    // Insert gaps into the alignment based on the cigar vector
    // -------------------------------------------------------------------------
    using std::get;
    auto current_ref_pos = std::ranges::begin(get<0>(alignment));
    auto current_read_pos = std::ranges::begin(get<1>(alignment));

    for (auto [cigar_count, cigar_operation] : cigar_vector)
    {
        if (('M'_cigar_operation == cigar_operation) || ('='_cigar_operation == cigar_operation)
            || ('X'_cigar_operation == cigar_operation))
        {
            std::ranges::advance(current_ref_pos, cigar_count);
            std::ranges::advance(current_read_pos, cigar_count);
        }
        else if (('D'_cigar_operation == cigar_operation) || ('N'_cigar_operation == cigar_operation))
        {
            // insert gaps into query
            current_read_pos = get<1>(alignment).insert_gap(current_read_pos, cigar_count);
            ++current_read_pos;
            std::ranges::advance(current_ref_pos, cigar_count);
        }
        else if (('I'_cigar_operation == cigar_operation)) // Insert gaps into ref
        {
            current_ref_pos = get<0>(alignment).insert_gap(current_ref_pos, cigar_count);
            ++current_ref_pos;
            std::ranges::advance(current_read_pos, cigar_count);
        }
        else if (('P'_cigar_operation == cigar_operation)) // skip padding
        {
            current_ref_pos = get<0>(alignment).insert_gap(current_ref_pos, cigar_count);
            ++current_ref_pos;

            current_read_pos = get<1>(alignment).insert_gap(current_read_pos, cigar_count);
            ++current_read_pos;
        }
        // S and H are ignored as they are handled by cropping the sequence
    }

    return alignment;
}

/*!\overload
 * \ingroup cigar_conversion
 */
template <typename reference_type, typename sequence_type>
inline auto alignment_from_cigar(std::string const & cigar_string,
                                 reference_type const & reference,
                                 uint32_t const zero_based_reference_start_position,
                                 sequence_type const & query)
{
    alignment_from_cigar(detail::parse_cigar(cigar_string), reference, zero_based_reference_start_position, query);
}

} // namespace seqan3
