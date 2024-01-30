// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::detail::all.
 */

#pragma once

#include <seqan3/contrib/std/all_view.hpp>
#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief A move-only view that takes unique ownership of a range.
 * \ingroup core
 * \sa https://en.cppreference.com/w/cpp/ranges/owning_view
 */
using SEQAN3_DOXYGEN_ONLY(owning_view =) seqan::stl::ranges::owning_view;

/*!\brief Returns a view that includes all elements of the range argument.
 * \ingroup core
 * \sa https://en.cppreference.com/w/cpp/ranges/all_view
 * \details
 * This implements the new C++20 behaviour that is only available with gcc12 and newer.
 * In contrast to the old std::views::all, rvalue ranges can be bound.
 * \returns
 *   * `rng` if `rng` is a view.
 *   * A std::ranges::ref_view of `rng` if that expression is valid.
 *   * Otherwise, a seqan3::detail::owning_view of `rng`.
 */
using SEQAN3_DOXYGEN_ONLY(all =) seqan::stl::views::all;

/*!\brief Returns the type that results from appying seqan3::detail::all to a range.
 * \ingroup core
 */
using SEQAN3_DOXYGEN_ONLY(all_t =) seqan::stl::views::all_t;

} // namespace seqan3::detail
