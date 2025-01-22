// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link io_sam_file IO / SAM File submodule \endlink.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

/*!\defgroup io_sam_file SAM File
 * \ingroup io
 * \brief Provides files and formats for handling read mapping data.
 *
 * # Introduction
 *
 * SAM/BAM files are primarily used to store pairwise alignments of read mapping data.
 *
 * \note For a step-by-step guide take a look at our tutorial: \ref tutorial_sam_file.
 *
 * \copydetails seqan3::sam_file_input::field_ids
 *
 * All of these fields are retrieved by default (and in that order).
 *
 * \include{doc} doc/fragments/io_sam_file_input.md
 *
 * \include{doc} doc/fragments/io_sam_file_output.md
 */

#pragma once

#include <seqan3/io/sam_file/format_bam.hpp>
#include <seqan3/io/sam_file/format_sam.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/input_format_concept.hpp>
#include <seqan3/io/sam_file/input_options.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/output_format_concept.hpp>
#include <seqan3/io/sam_file/output_options.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
