// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::sam_file_header class.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <deque>
#include <ranges>
#include <unordered_map>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

namespace seqan3
{

/*!\brief Stores information of the program/tool that was used to create a SAM/BAM file.
 * \ingroup io_sam_file
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
struct sam_file_program_info_t
{
    std::string id;                //!< A unique (file scope) id.
    std::string name;              //!< The official name.
    std::string command_line_call; //!< The command line call that produces the file.
    std::string previous;          //!< The id of the previous program if program calls were chained.
    std::string description;       //!< A description of the program and/or program call.
    std::string version;           //!< The program/tool version.
    std::string user_tags;         //!< Additional user-defined tags.
};

/*!\brief Stores the header information of SAM/BAM files.
 * \ingroup io_sam_file
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
template <std::ranges::forward_range ref_ids_type = std::deque<std::string>>
class sam_file_header
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    sam_file_header() = default;                                    //!< Defaulted.
    sam_file_header(sam_file_header const &) = default;             //!< Defaulted.
    sam_file_header & operator=(sam_file_header const &) = default; //!< Defaulted.
    sam_file_header(sam_file_header &&) = default;                  //!< Defaulted.
    sam_file_header & operator=(sam_file_header &&) = default;      //!< Defaulted.
    ~sam_file_header() = default;                                   //!< Defaulted.

    /*!\brief Construct from a range of reference ids.
     * \param[in] ref_ids The range over reference ids.
     */
    sam_file_header(ref_ids_type ref_ids) : reference_ids{std::move(ref_ids)}
    {}
    //!\}

    using program_info_t =
        sam_file_program_info_t; //!< Stores information of the program/tool that was used to create the file.
    std::string format_version;  //!< The file format version. Note: this is overwritten by our formats on output.
    std::string sorting;         //!< The sorting of the file. SAM: [unknown, unsorted, queryname, coordinate].
    std::string
        subsorting; //!< The sub-sorting of the file. SAM: [unknown, unsorted, queryname, coordinate](:[A-Za-z0-9_-]+)+.
    std::string grouping; //!< The grouping of the file. SAM: [none, query, reference].

    std::vector<program_info_t> program_infos; //!< The list of program information.

    std::vector<std::string> comments; //!< The list of comments.

private:
    //!\brief The key's type of ref_dict.
    using key_type = std::conditional_t<std::ranges::contiguous_range<std::ranges::range_reference_t<ref_ids_type>>,
                                        std::span<range_innermost_value_t<ref_ids_type> const>,
                                        type_reduce_t<std::ranges::range_reference_t<ref_ids_type>>>;
    //!\brief The reference ids.
    ref_ids_type reference_ids{};

    //!\brief Custom hash function since std::hash is not defined for all range types (e.g. std::span<char>).
    struct key_hasher
    {
        //!\brief Hash a key.
        template <typename key_t>
        size_t operator()(key_t && key) const noexcept
        {
            using char_t = std::ranges::range_value_t<key_t>;
            size_t result{0};
            std::hash<char_t> h{};
            for (char_t character : key)
            {
                result *= 0x8F'3F'73'B5'CF'1C'9A'DE;
                result += h(character);
            }
            return result;
        }
    };

public:
    /*!\brief The range of reference ids.
     *
     * \details
     *
     * This member function gives you access to the range of reference ids.
     *
     * When reading a file, there are three scenarios:
     * 1) Reference id information is provided on construction. In this case, no copy is made but this function
     *    gives you a reference to the provided range. When reading the header or the records, their reference
     *    information will be checked against the given input.
     * 2) No reference information is provided on construction but the \@SQ tags are present in the header.
     *    In this case, the reference id information is extracted from the header and this member function provides
     *    access to them. When reading the records, their reference id information will be checked against the header
     *    information.
     * 3) No reference information is provided on construction an no \@SQ tags are present in the header.
     *    In this case, the reference information is parsed from the records field::ref_id and stored in the header.
     *    This member function then provides access to the unique list of reference ids encountered in the records.
     */
    ref_ids_type & ref_ids()
    {
        return reference_ids;
    }

    /*!\brief The reference information. (used by the SAM/BAM format)
     *
     * \details
     *
     * The reference information store the length (\@LN tag) and
     * additional information of each reference sequence in the file. The record
     * must then store only the index of the reference.
     * The name and length information are required if the header is provided
     * and each reference sequence that is referred to in any of the records
     * must be present in the dictionary, otherwise a seqan3::format_error will
     * be thrown upon reading or writing a file.
     *
     * The additional information (2nd tuple entry) must model
     * the following formatting rules: The information is given in a tab separated
     * TAG:VALUE format, where TAG must be one of [AH, AN, AS, m5, SP, UR].
     * The following information and rules apply for each tag (taken from the SAM specs):
     *
     * * **AH:** Indicates that this sequence is an alternate locus. The value is the locus in the primary assembly for
     *           which this sequence is an alternative, in the format 'chr:start-end', 'chr' (if known), or '*' (if
     *           unknown), where 'chr' is a sequence in the primary assembly. Must not be present on sequences in the
     *           primary assembly.
     * * **AN:** Alternative reference sequence names. A comma-separated list of alternative names that tools may use
     *           when referring to this reference sequence. These alternative names are not used elsewhere within the
     *           SAM file; in  particular, they must not appear in SAM records’ RNAME or RNEXT fields. regular
     *           expression : name (, name )* where name is [0-9A-Za-z][0-9A-Za-z*+.@ \|-]*
     * * **AS:** Genome assembly identifier.
     * * **M5:** MD5 checksum of the sequence.  See Section 1.3.1
     * * **SP:** Species.
     * * **UR:** URI of the sequence.  This value may start with one of the standard protocols, e.g http:  or ftp:. If
     *           it does not start with one of these protocols, it is assumed to be a file-system path
     */
    std::vector<std::tuple<int32_t, std::string>> ref_id_info{};

    //!\brief The mapping of reference id to position in the ref_ids() range and the ref_id_info range.
    std::unordered_map<key_type, int32_t, key_hasher, detail::view_equality_fn> ref_dict{};

    /*!\brief The Read Group Dictionary (used by the SAM/BAM format).
     *
     * \details
     *
     * The read group dictionary stores the group id and
     * additional information of each read group in the file. The record
     * may store a RG tag information referencing one of the stored id's.
     * The id information is required if the header is provided.
     *
     * The additional information (2nd tuple entry) for the SAM format must follow
     * the following formatting rules: The information is given in a tab separated
     * TAG:VALUE format, where TAG must be one of [AH, AN, AS, m5, SP, UR].
     * The following information and rules apply for each tag (taken from the SAM specs):
     *
     * * **BC:** Barcode sequence identifying the sample or library. This value is the expected barcode bases as read by
     *           the sequencing machine in the absence of errors. If there are several barcodes for the sample/library
     *           (e.g., one on each end of the template), the recommended implementation concatenates all the barcodes
     *           separating them with hyphens ('-').
     * * **CN:** Name of sequencing center producing the read.
     * * **DS:** Description.  UTF-8 encoding may be used.
     * * **DT:** Date the run was produced (ISO8601 date or date/time).
     * * **FO:** Flow order. The array of nucleotide bases that correspond to the nucleotides used for each flow of each
     *           read. Multi-base flows are encoded in IUPAC format, and non-nucleotide flows by various other
     *           characters. Format : /\*\|[ACMGRSVTWYHKDBN]+/
     * * **KS:** The array of nucleotide bases that correspond to the key sequence of each read.
     * * **LB:** Library.
     * * **PG:** Programs used for processing the read group.
     * * **PI:** Predicted median insert size.
     * * **PL:** Platform/technology used to produce the reads. Valid values : CAPILLARY, LS454, ILLUMINA, SOLID,
     *           HELICOS, IONTORRENT, ONT, and PACBIO.
     * * **PM:** Platform model. Free-form text providing further details of the platform/technology used.
     * * **PU:** Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.
     * * **SM:** Sample. Use pool name where a pool is being sequenced.
     */
    std::vector<std::pair<std::string, std::string>> read_groups;

    std::string user_tags; //!< Additional user-defined tags.
};

} // namespace seqan3
