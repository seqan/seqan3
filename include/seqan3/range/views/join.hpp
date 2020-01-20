// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::join.
 */

#pragma once

#include <range/v3/view/join.hpp>

#include <seqan3/core/platform.hpp>
#include <seqan3/std/ranges>

namespace seqan3::views
{

//TODO reimplement me and rename to join_with
using ::ranges::views::join;

} // namespace seqan3::views
