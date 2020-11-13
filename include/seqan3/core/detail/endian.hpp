// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::detail::endian.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \deprecated This header is deprecated and will be removed in SeqAn-3.1.0;. Please use std::endian and \#include
 *             <seqan3/std/bit> instead.
 */

#pragma once

#include <seqan3/std/bit>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

    //!\brief Alias for std::endian.
    using endian [[SEQAN3_DEPRECATED_310]] = std::endian;

} // namespace seqan3::detail

SEQAN3_DEPRECATED_HEADER(
   "This header is deprecated and will be removed in SeqAn-3.1.0; Please use std::endian and #include <seqan3/std/bit> instead.")
