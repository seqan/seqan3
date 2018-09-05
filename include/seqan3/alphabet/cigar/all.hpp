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
 * \brief Meta-header for the cigar submodule.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>

/*!\defgroup cigar Cigar
 * \ingroup alphabet
 * \brief Contains the cigar_op alphabet, cigar alphabet and the cigar_vector.
 *
 * ### Introduction
 * The cigar string is a way to represent an alignment of two
 * sequences and was introduced by the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).
 * It is most commonly used for representing the alignment between a read and
 * the reference.
 *
 * A cigar string may look like this: `20M2D20M1I10M`, which means that the first
 * 20 bases are aligned (either by match or mismatch), then 2 bases are deleted,
 * again 20 bases are aligned, then 1 base is inserted, and finally 10 aligned.
 * This shows that a cigar string consist of a series of cigar elements, that
 * have an operator value (e.g. 'M') and a length value (e.g. 20).
 *
 * The possible operator values, modelled in the seqan3::cigar_op alphabet are
 * the following:
 *
 * | Operator | Description |
 * |----------|-------------|
 * | M        | Match/Aligned Bases; can be either an alignment match or mismatch of two aligned bases. |
 * | D        | Deletion; the nucleotide is present in the reference but not in the read. |
 * | I        | Insertion; the nucleotide is present in the read  but not in the reference. |
 * | S        | Soft Clipping; the clipped nucleotides are present in the read sequence but not part of the alignment. |
 * | H        | Hard Clipping; the clipped nucleotides are not present in the read. |
 * | N        | Skipped region; a region of nucleotides is not present in the read |
 * | P        | Padding; padded area in the read and not in the reference |
 * | X        | Alignment Mismatch; the aligned characters are not equal. |
 * | =        | Alignment Match; the aligned characters are equal. |
 *
 * A cigar element is modelled by the seqan3::cigar alphabet combining a
 * seqan3::cigar_op character with a length value.
 *
 * Some example code for the usage of seqan3::cigar_op and seqan3::cigar:
 *
 * \snippet test/snippet/alphabet/cigar/cigar_op.cpp general
 */
