// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-include for the alignment IO submodule.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

/*!\defgroup alignment_file Alignment File
 * \ingroup io
 * \brief Provides files and formats for handling alignment data.
 *
 * ### Introduction
 *
 * Alignment files are primarily used to store pairwise alignments of two biological sequences and often come with
 * many additional information. Well-known formats include the SAM/BAM format used to store read mapping data or the
 * BLAST format that stores the results of a query search against a data base.
 *
 * \note For a step-by-step guide take a look at our tutorial: \ref tutorial_alignment_file.
 *
 * \copydetails seqan3::alignment_file_input::field_ids
 *
 * All of these fields are retrieved by default (and in that order).
 * Note that some of the fields are specific to the SAM format (e.g. seqan3::field::flag) while others are specific to
 * BLAST format (e.g. seqan3::field::bit_score). Please see the corresponding formats for more details.
 */

#include <seqan3/io/alignment_file/format_bam.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/input_options.hpp>
#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
