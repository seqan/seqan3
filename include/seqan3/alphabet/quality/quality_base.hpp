// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::phred42 quality scores.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/alphabet/quality/phred_base.hpp>
 *             instead.
 */

#pragma once

#include <seqan3/alphabet/quality/phred_base.hpp>

namespace seqan3
{
//!\deprecated Please use seqan3::phred_base instead.
template <typename derived_type, size_t size>
using quality_base SEQAN3_DEPRECATED_310 = seqan3::phred_base<derived_type, size>;
} // namespace seqan3

SEQAN3_DEPRECATED_HEADER(
   "This header is deprecated and will be removed in SeqAn-3.1.0; Please #include <seqan3/alphabet/quality/phred_base.hpp> instead.")
