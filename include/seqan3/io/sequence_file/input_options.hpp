// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sequence_file_input_options.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \ingroup io_sequence_file
 * \tparam sequence_legal_alphabet_ The sequence legal alphabet exposed as type trait to the format.
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */
template <typename sequence_legal_alphabet>
struct sequence_file_input_options
{
    //!\brief Read the ID string only up until the first whitespace character.
    bool truncate_ids = false;
    //!\brief Read the complete_header into the seqan3::field::id for embl or genbank format.
    bool embl_genbank_complete_header = false;
    //!\brief Remove spaces after ">" (or ";") before the actual ID.
    bool fasta_ignore_blanks_before_id = true;
};

} // namespace seqan3
