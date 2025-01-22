// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sequence_file_output_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \ingroup io_sequence_file
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */
struct sequence_file_output_options
{
    //!\brief Begin the ID line with ";" instead of ">" (not recommended).
    bool fasta_legacy_id_marker = false;
    //!\brief Insert a single space after ">" (or ";") before the actual ID.
    bool fasta_blank_before_id = false;
    //!\brief Inserts linebreaks after every n-th letter in the sequence; 0 means no linebreaks.
    uint32_t fasta_letters_per_line = 80;
    //TODO:
    //     bool        fasta_charcounts        = false;

    //!\brief Whether to write the ID only '@' or also after '+' line.
    bool fastq_double_id = false;

    /*!\brief The default plain text line-ending is "\n", but on Windows an additional carriage return is
     *        recommended ("\r\n" for line-ending).
     */
    bool add_carriage_return = false;

    //!\brief Complete header given for embl or genbank
    bool embl_genbank_complete_header = false;
};

} // namespace seqan3
