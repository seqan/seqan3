// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Meta-header for the cigar submodule; includes all headers from alphabet/cigar/.
 */

 #pragma once

 #include <seqan3/alphabet/cigar/cigar_op.hpp>

 /*!\defgroup cigar CIGAR
  * \brief Provides (semi-)alphabets for representing elements in CIGAR strings.
  * \ingroup alphabet
  *
  * \details
  *
  * ###Introduction
  *
  * CIGAR strings are combinations of count values and CIGAR operations, representing an alignment as a sequence of
  * edit operations. This submodule has two different
  * alphabets. One is the seqan3::cigar_op alphabet, which is a base seqan3::alphabet implementation. This
  * contains all valid symbols contained in CIGAR strings. The other alphabet is the seqan3::cigar alphabet, which
  * is an alphabet tuple. It combines the seqan3::cigar_op alphabet with a count value,
  * such that one can represent an entire CIGAR string with a std::vector of seqan3::cigar values.
  *
  * The following table outlines the valid characters in the seqan3::cigar_op alphabet.
  *
  * | Letter | Description                                                                                     |
  * |--------|-------------------------------------------------------------------------------------------------|
  * | M      | Alignment match (can be a sequence match or mismatch, used only in basic CIGAR representations) |
  * | I      | Insertion to the reference                                                                      |
  * | D      | Deletion from the reference                                                                     |
  * | N      | Skipped region from the reference                                                               |
  * | S      | Soft clipping (clipped sequences present in SEQ)                                                |
  * | H      | Hard clipping (clipped sequences NOT present in SEQ)                                            |
  * | P      | Padding (silent deletion from padded reference)                                                 |
  * | =      | Sequence match                                                                                  |
  * | X      | Sequence mismatch                                                                               |
  */
