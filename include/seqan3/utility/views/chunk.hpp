// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::views::chunk.
 */

#pragma once

#include <seqan3/contrib/std/chunk_view.hpp>
#include <seqan3/core/platform.hpp>

namespace seqan3::views
{

/*!\brief A view adaptor that divides a range into chunks.
 * \ingroup utility_views
 * \noapi{This is a implementation of the C++23 chunk_view. It will be replaced with std::views::chunk.}
 * \sa https://en.cppreference.com/w/cpp/ranges/chunk_view
 */
using SEQAN3_DOXYGEN_ONLY(chunk =) seqan::stl::views::chunk;

} // namespace seqan3::views
