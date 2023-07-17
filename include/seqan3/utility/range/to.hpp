// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Simon Gene Gottlieb <simon.gottlieb AT fu-berlin.de>
 * \brief Provides seqan3::ranges::to.
 */

#pragma once

#include <seqan3/contrib/std/to.hpp>
#include <seqan3/core/platform.hpp>

namespace seqan3::ranges
{

/*!\brief Converts a range to a container.
 * \ingroup utility_range
 * \noapi{This is a implementation of the C++23 ranges::to. It will be replaced with std::ranges::to.}
 * \sa https://en.cppreference.com/w/cpp/ranges/to
 */
using SEQAN3_DOXYGEN_ONLY(to =) seqan::std::ranges::to;

} // namespace seqan3::ranges
