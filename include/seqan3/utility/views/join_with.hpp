// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::views::join_with.
 */

#pragma once

#include <seqan3/contrib/std/join_with_view.hpp>
#include <seqan3/core/platform.hpp>

namespace seqan3::views
{

/*!\brief A view adaptor that represents view consisting of the sequence obtained from flattening a view of ranges, with
 *        every element of the delimiter inserted in between elements of the view. The delimiter can be a single element
 *        or a view of elements.
 * \ingroup utility_views
 * \noapi{This is a implementation of the C++23 join_with_view. It will be replaced with ::std::views::join_with.}
 * \sa https://en.cppreference.com/w/cpp/ranges/join_with_view
 */
using SEQAN3_DOXYGEN_ONLY(join_with =) seqan::stl::views::join_with;

} // namespace seqan3::views
