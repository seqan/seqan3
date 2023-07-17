// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Meta-header for the \link alphabet_cigar Alphabet / CIGAR submodule \endlink.
 */

/*!\defgroup alphabet_cigar CIGAR
  * \brief Provides (semi-)alphabets for representing elements in CIGAR strings.
  * \ingroup alphabet
  * \see alphabet
  *
  * \details
  *
  * ###Introduction
  *
  * CIGAR strings are combinations of count values and CIGAR operations, representing an alignment as a sequence of
  * edit operations. This submodule has two different
  * alphabets. One is the seqan3::cigar::operation alphabet, which is a base seqan3::alphabet implementation. This
  * contains all valid symbols contained in CIGAR strings. The other alphabet is the seqan3::cigar alphabet, which
  * is an alphabet tuple. It combines the seqan3::cigar::operation alphabet with a count value,
  * such that one can represent an entire CIGAR string with a std::vector of seqan3::cigar values.
  *
  * The following table outlines the valid characters in the seqan3::cigar::operation alphabet.
  *
  * \include{doc} doc/fragments/cigar_operation_table.md
  */

#pragma once

#include <seqan3/alphabet/cigar/cigar.hpp>
