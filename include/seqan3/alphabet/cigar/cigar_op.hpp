// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Introduces the seqan3::cigar_op alphabet.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/alphabet/cigar/cigar.hpp>
 *             instead.
 */

#pragma once

#include <seqan3/alphabet/cigar/cigar.hpp>

namespace seqan3
{
//!\deprecated Please use seqan3::cigar::operation instead.
using cigar_op SEQAN3_DEPRECATED_310 = seqan3::cigar::operation;
} // namespace seqan3

SEQAN3_DEPRECATED_HEADER(
   "This header is deprecated and will be removed in SeqAn-3.1.0; Please #include <seqan3/alphabet/cigar/cigar.hpp> instead.")
