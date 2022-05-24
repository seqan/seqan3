// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides seqan3::views::chunk.
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/core/platform.hpp>

#include <range/v3/view/chunk.hpp>

namespace seqan3::views
{

/*!\brief A chunk view
 * \ingroup utility_views
 * \details
 * \noapi{This is currently range-v3's chunk implementation.}
 */
inline constexpr auto chunk = ::ranges::views::chunk;

} // namespace seqan3::views
