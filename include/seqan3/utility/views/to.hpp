// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::to.
 */

#pragma once

#include <seqan3/std/ranges>

#include <range/v3/range/conversion.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::views
{

/*!\brief A to view
 * \ingroup views
 * \details
 * \noapi{This is currently range-v3's to implementation.}
 */
#if SEQAN3_DOXYGEN_ONLY(1)0
inline constexpr auto to;
#else // ^^^ doxygen only / real import vvv
using ::ranges::to;
#endif // SEQAN3_DOXYGEN_ONLY(1)0

} // namespace seqan3::views
