// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides helper data structures for the seqan3::sam_file_output.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>

namespace seqan3
{

//!\brief Type tag which indicates that no reference information has been passed to the SAM file on construction.
//!\ingroup io_sam_file
struct ref_info_not_given
{};

/*!\brief An enum flag that describes the properties of an aligned read (given as a SAM record).
 * \ingroup io_sam_file
 * \implements seqan3::enum_bitwise_operators
 * \sa seqan3::enum_bitwise_operators enables combining enum values.
 *
 * The SAM flag are bitwise flags, which means that each value corresponds to a specific bit that is set and that they
 * can be combined and tested using binary operations.
 * See this [tutorial](https://www.codeproject.com/Articles/13740/The-Beginner-s-Guide-to-Using-Enum-Flags) for an
 * introduction on bitwise operations on enum flags.
 *
 * Example:
 *
 * \include test/snippet/io/sam_file/sam_flags.cpp
 *
 * Adapted from the [SAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) are the following additional
 * information to some flag values:
 * * For each read/contig in a SAM file, it is required that one and only one line associated with the read
 *   has neither the seqan3::sam_flag::secondary_alignment nor the seqan3::sam_flag::supplementary_alignment flag value
 *   set (satisfies `FLAG & 0x900 == 0 `). This line is called the **primary alignment** of the read.
 * * seqan3::sam_flag::secondary_alignment (bit `0x100`) marks the alignment not to be used in certain analyses when
 *   the tools in use are aware of this bit. It is typically used to flag alternative mappings when multiple mappings
 *   are presented in a SAM.
 * * seqan3::sam_flag::supplementary_alignment (bit `0x800`) indicates that the corresponding alignment line is part
 *   of a chimeric alignment. If the SAM/BAM file corresponds to long reads (nanopore/pacbio) this happens when
 *   reads are split before being aligned and the best matching part is marked as primary, while all other aligned
 *   parts are marked supplementary.
 * * seqan3::sam_flag::unmapped (bit `0x4`) is the only reliable place to tell whether the read is unmapped.
 *   If seqan3::sam_flag::unmapped is set, no assumptions can be made about RNAME, POS, CIGAR, MAPQ, and
 *   seqan3::sam_flag::proper_pair, seqan3::sam_flag::secondary_alignment, and seqan3::sam_flag::supplementary_alignment
 *   (bits `0x2`, `0x100`, and `0x800`).
 * * seqan3::sam_flag::on_reverse_strand (bit `0x10`) indicates  whether the read sequence has been reverse complemented
 *   and the quality string is reversed.  When bit seqan3::sam_flag::unmapped (`0x4`) is unset, this
 *   corresponds to the strand to which the segment has been mapped: seqan3::sam_flag::on_reverse_strand (bit `0x10`)
 *   unset indicates the forward strand, while set indicates the reverse strand. When seqan3::sam_flag::unmapped (`0x4`)
 *   is set, this indicates whether the unmapped read is stored in its original orientation as it came off the
 *   sequencing machine.
 * * seqan3::sam_flag::first_in_pair and seqan3::sam_flag::second_in_pair (bits `0x40` and `0x80`) reflect the read
 *   ordering within each template inherent in the sequencing technology used. If seqan3::sam_flag::first_in_pair and
 *   seqan3::sam_flag::second_in_pair (`0x40` and `0x80`) are both set, the read is part of a linear template, but it
 *   is neither the first nor the last read. If both are unset, the index of the read in the template is unknown.
 *   This may happen for a non-linear template or when this information is lost during data processing.
 * * If seqan3::sam_flag::paired (bit `0x1`) is unset, no assumptions can be made about seqan3::sam_flag::proper_pair,
 *   seqan3::sam_flag::mate_unmapped, seqan3::sam_flag::mate_on_reverse_strand, seqan3::sam_flag::first_in_pair and
 *   seqan3::sam_flag::second_in_pair (bits `0x2`, `0x8`, `0x20`, `0x40` and `0x80`).
 *
 * \sa https://broadinstitute.github.io/picard/explain-flags.html
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
enum class sam_flag : uint16_t
{
    none = 0,                       //!< None of the flags below are set.
    paired = 0x1,                   //!< The aligned read is paired (paired-end sequencing).
    proper_pair = 0x2,              //!< The two aligned reads in a pair have a proper distance between each other.
    unmapped = 0x4,                 //!< The read is not mapped to a reference (unaligned).
    mate_unmapped = 0x8,            //!< The mate of this read is not mapped to a reference (unaligned).
    on_reverse_strand = 0x10,       //!< The read sequence has been reverse complemented before being mapped (aligned).
    mate_on_reverse_strand = 0x20,  //!< The mate sequence has been reverse complemented before being mapped (aligned).
    first_in_pair = 0x40,           //!< Indicates the ordering (see details in the seqan3::sam_flag description).
    second_in_pair = 0x80,          //!< Indicates the ordering (see details in the seqan3::sam_flag description).
    secondary_alignment = 0x100,    //!< This read alignment is an alternative (possibly suboptimal) to the primary.
    failed_filter = 0x200,          //!< The read alignment failed a filter, e.g. quality controls.
    duplicate = 0x400,              //!< The read is marked as a PCR duplicate or optical duplicate.
    supplementary_alignment = 0x800 //!< This sequence is part of a split alignment and is not the primary alignment.
};

//!\cond DEV
//!\brief Enables bitwise operations for seqan3::sam_flags.
//!\ingroup io_sam_file
//!\sa seqan3::enum_bitwise_operators enables combining enum values.
template <>
inline constexpr bool add_enum_bitwise_operators<sam_flag> = true;
//!\endcond

/*!\brief A sam_flag can be printed as an integer value.
 * \ingroup io_sam_file
 */
template <>
struct sam_flag_printer<sam_flag>
{
    /*!\brief Prints the sam flag.
     * \tparam stream_t The output stream type.
     * \param[in,out] stream The output stream.
     * \param[in] arg The sam flag to print.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, sam_flag const arg) const
    {
        stream << static_cast<int16_t>(arg);
    }
};

} // namespace seqan3
