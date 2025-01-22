// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link alignment Alignment module \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup alignment Sequence Alignment
 * \brief The alignment module contains concepts, algorithms and classes that are related to the computation of
 *        pairwise and multiple sequence alignments.
 *
 * \details
 *
 * # Sequence Alignment
 *
 * In bioinformatics, a sequence alignment is a way of arranging the sequences of DNA, RNA, or protein to identify
 * regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between
 * the sequences. Aligned sequences of nucleotide or amino acid residues are typically represented as rows within
 * a matrix. Gaps are inserted between the residues so that identical or similar characters are aligned in successive
 * columns. Sequence alignments are also used for non-biological sequences, such as calculating the distance cost
 * between strings in a natural language or in financial data. [1]
 *
 * # Pairwise Sequence Alignment
 *
 * SeqAn offers a generic multi-purpose alignment library comprising all widely known alignment algorithms as well as
 * many special algorithms. These algorithms are all accessible through an easy to use alignment interface which
 * is described in \ref alignment_pairwise.
 *
 * The following code snippet demonstrates a simple use of the pairwise alignment interface.
 *
 * \include doc/tutorial/08_pairwise_alignment/pairwise_alignment_first_global.cpp
 *
 * # Multiple Sequence Alignment
 *
 * The current version of SeqAn does not offer multiple sequence alignments (MSA). Please reach out to us with
 * a specific use case we should consider in future versions.
 *
 * # Alignments represented as CIGAR String used in SAM/BAM Files
 *
 * A common file format to store (semi) alignments is the SAM/BAM format. In a SAM/BAM file, the alignment is
 * represented as a CIGAR string. To allow back and forth conversion from a CIGAR string to the alignment
 * representation in SeqAn, we provide the following functions:
 *
 * * seqan3::alignment_from_cigar
 * * seqan3::cigar_from_alignment
 *
 * For reading and writing SAM/BAM files, we provide the seqan3::sam_file_input and seqan3::sam_file_ouput.
 *
 * # References
 *
 * [1] https://en.wikipedia.org/wiki/Sequence_alignment
 */

#pragma once

#include <seqan3/alignment/aligned_sequence/all.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/decorator/all.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/matrix/all.hpp>
#include <seqan3/alignment/pairwise/all.hpp>
#include <seqan3/alignment/scoring/all.hpp>
