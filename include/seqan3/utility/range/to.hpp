// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
using SEQAN3_DOXYGEN_ONLY(to =) seqan::stl::ranges::to;

} // namespace seqan3::ranges
