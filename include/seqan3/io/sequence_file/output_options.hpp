// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::sequence_file_output_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

//!\brief The options type defines various option members that influence the behaviour of all or some formats.
struct sequence_file_output_options
{
    //!\brief Begin the ID line with ";" instead of ">" (not recommended).
    bool        fasta_legacy_id_marker  = false;
    //!\brief Insert a single space after ">" (or ";") before the actual ID.
    bool        fasta_blank_before_id   = true;
    //!\brief Inserts linebreaks after every n-th letter in the sequence; 0 means no linebreaks.
    uint32_t    fasta_letters_per_line  = 80;
    //TODO:
//     bool        fasta_charcounts        = false;


    //!\brief Whether to write the ID only '@' or also after '+' line.
    bool        fastq_double_id         = false;

    /*!\brief The default plain text line-ending is "\n", but on Windows an additional carriage return is
     *        recommended ("\r\n" for line-ending).
     */
    bool        add_carriage_return     = false;

    //!\brief Complete header given for embl or genbank
    bool        embl_genbank_complete_header  = false;
};

} // namespace seqan3
