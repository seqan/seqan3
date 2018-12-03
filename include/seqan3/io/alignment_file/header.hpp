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
 * \brief Provides the seqan3::alignment_file_header class.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <unordered_map>
#include <vector>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

//!\brief Stores the header information of alignment files.
//!\ingroup alignment_file
struct alignment_file_header
{
    //!\brief Stores information of the program/tool that was used to create the file.
    struct program_info_t
    {
        std::string id;                //!< A unique (file scope) id.
        std::string name;              //!< The official name.
        std::string command_line_call; //!< The command line call that produces the file.
        std::string previous;          //!< The id of the previous program if program calls were chained.
        std::string description;       //!< A description of the program and/or program call.
        std::string version;           //!< The program/tool version.
    };

    std::string format_version;     //!< The file format version. Note: this is overwritten by our formats on output.
    std::string sorting{"unknown"}; //!< The sorting state of the file. SAM: [unknown, unsorted, queryname, coordinate].
    std::string grouping{"none"};   //!< The grouping state of the file. SAM: [none, query, reference].

    std::vector<program_info_t> program_infos; //!< The list of program information.

    std::vector<std::string> comments;         //!< The list of comments.

    /*!\brief The Reference Dictionary. (used by the SAM/BAM format)
     *
     * \details
     *
     * The reference dictionary stores the reference name, its length and
     * additional information of each reference sequence in the file. The record
     * may (SAM) or must (BAM) then store only the index of the reference.
     * The name and length information are required if the header is provided
     * and each reference sequence that is referred to in any of the records
     * must be present in the dictionary, otherwise a seqan3::format_error will
     * be thrown upon reading or writing a file.
     *
     * The additional information (2nd tuple entry) for the SAM format must follow
     * the following formatting rules: The information is given in a tab separated
     * TAG:VALUE format, where TAG must be one of [AH, AN, AS, m5, SP, UR].
     * The following information and rules apply for each tag (taken from the SAM specs):
     *
     * | TAG | Description and Rules                                            |
     * | --- | ---------------------------------------------------------------- |
     * | AH  | Indicates that this sequence is an alternate locus. The value is
     *         the locus in the primary assembly for which this sequence is an
     *         alternative, in the format 'chr:start-end', 'chr' (if known),
     *         or '*' (if unknown), where 'chr' is a sequence in the primary assembly.
     *         Must not be present on sequences in the primary assembly. |
     * | AN  | Alternative reference sequence names. A comma-separated list
     *         of alternative names that tools may use when referring to this
     *         reference sequence. These alternative names are not used elsewhere
     *         within the SAM file; in  particular, they must not appear in alignment
     *         records’ RNAME or RNEXT fields. Regular expression : name (, name )*
     *         where name is [0-9A-Za-z][0-9A-Za-z*+.@ |-]* |
     * | AS  | Genome assembly identifier. |
     * | M5  | MD5 checksum of the sequence.  See Section 1.3.1 |
     * | SP  | Species. |
     * | UR  | URI of the sequence.  This value may start with one of the standard
     *         protocols, e.g http:  or ftp:. If it does not start with one of these
     *         protocols, it is assumed to be a file-system path |
     */
    std::unordered_map<std::string, std::tuple<uint32_t, std::string>> ref_dict;

    /*!\brief The Read Group Dictionary. (used by the SAM/BAM format)
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
     * | TAG | Description and Rules                                            |
     * | --- | ---------------------------------------------------------------- |
     * | BC  | Barcode sequence identifying the sample or library. This value is
     *         the expected barcode bases as read by the sequencing machine in
     *         the absence of errors. If there are several barcodes for the
     *         sample/library (e.g., one on each end of the template), the
     *         recommended implementation concatenates all the barcodes separating
     *         them with hyphens ('-'). |
     * | CN  | Name of sequencing center producing the read. |
     * | DS  | Description.  UTF-8 encoding may be used. |
     * | DT  | Date the run was produced (ISO8601 date or date/time). |
     * | FO  | Flow order. The array of nucleotide bases that correspond to the
     *         nucleotides used for each flow of each read. Multi-base flows are
     *         encoded in IUPAC format, and non-nucleotide flows by various other
     *         characters. Format : /\*|[ACMGRSVTWYHKDBN]+/ |
     * | KS  | The array of nucleotide bases that correspond to the key sequence of each read. |
     * | LB  | Library. |
     * | PG  | Programs used for processing the read group. |
     * | PI  | Predicted median insert size. |
     * | PL  | Platform/technology used to produce the reads.
     *         Valid values : CAPILLARY, LS454, ILLUMINA, SOLID, HELICOS, IONTORRENT, ONT, and PACBIO. |
     * | PM  | Platform model. Free-form text providing further details of the platform/technology used. |
     * | PU  | Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier. |
     * | SM  | Sample. Use pool name where a pool is being sequenced. |
     */
    std::vector<std::pair<std::string, std::string>> read_groups;
};

} // namespace seqan3
