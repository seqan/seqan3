// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link io_sequence_file IO / Sequence File submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

/*!\defgroup io_sequence_file Sequence File
 * \brief Provides files and formats for handling sequence data.
 * \ingroup io
 *
 * \include{doc} doc/fragments/sequence_file_input.md
 *
 * \include{doc} doc/fragments/sequence_file_output.md
 *
 * \see io
 * \see \ref tutorial_sequence_file
 */

#pragma once

#include <seqan3/io/sequence_file/format_embl.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/format_genbank.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/sequence_file/record.hpp>
