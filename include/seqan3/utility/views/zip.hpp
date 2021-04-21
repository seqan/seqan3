// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::zip.
 */

#pragma once

#include <seqan3/std/ranges>

#include <range/v3/view/zip.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::views
{

/*!\brief A zip view
 * \ingroup views
 * \details
 * \noapi{This is currently range-v3's zip implementation.}
 */
inline constexpr auto zip = ::ranges::views::zip;

} // namespace seqan3::views
