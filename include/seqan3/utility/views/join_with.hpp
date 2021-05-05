// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::join_with.
 */

#pragma once

#include <seqan3/std/ranges>

#include <range/v3/view/join.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::views
{

/*!\brief A join view.
 * \ingroup views
 * \deprecated Please use std::views::join or seqan3::views::join_with (if a separator is needed)
 */
SEQAN3_DEPRECATED_310 inline constexpr auto join = ::ranges::views::join;

/*!\brief A join view, please use std::views::join if you don't need a separator.
 * \ingroup views
 * \details
 * \noapi{This is currently range-v3's join implementation.}
 */
inline constexpr auto join_with = ::ranges::views::join;

} // namespace seqan3::views
