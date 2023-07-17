// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation helper for seqan3::views::zip and seqan3::views::join_with.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/core/platform.hpp>

//!\cond
namespace seqan3::detail::view_helper
{

template <class range_t>
concept simple_view = std::ranges::view<range_t> && std::ranges::range<range_t const>
                   && std::same_as<std::ranges::iterator_t<range_t>, std::ranges::iterator_t<range_t const>>
                   && std::same_as<std::ranges::sentinel_t<range_t>, std::ranges::sentinel_t<range_t const>>;

template <bool is_const, typename t>
using maybe_const = std::conditional_t<is_const, t const, t>;

} // namespace seqan3::detail::view_helper
//!\endcond
