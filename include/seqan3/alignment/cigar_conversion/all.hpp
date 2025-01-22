// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link cigar_conversion Alignment / CIGAR Conversion submodule \endlink.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

/*!\defgroup cigar_conversion CIGAR Conversion
 * \brief The CIGAR Conversion submodule contains utility functions to convert a CIGAR to an alignment or vice versa.
 * \ingroup alignment
 *
 * # Quick background on the CIGAR string
 *
 * The CIGAR string is a compact representation of an aligned read against a reference and was introduced by
 * the [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format. The SAM format stores the result of mapping
 * short/long read sequences from a sequencing experiment (e.g., Illumina/Nanopore) against a reference (e.g., hg38).
 *
 * # Creating an alignment from a CIGAR String
 *
 * To create an alignment from a CIGAR String, you can do the following:
 *
 * \include test/snippet/alignment/cigar_conversion/alignment_from_cigar_io.cpp
 *
 * # Creating a CIGAR string from an alignment
 *
 * To create a CIGAR String from an alignment, you can do the following:
 *
 * \include test/snippet/alignment/cigar_conversion/cigar_from_alignment_with_clipping.cpp
 *
 * \sa seqan3::alignment_from_cigar
 * \sa seqan3::cigar_from_alignment
 */

#pragma once

#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
