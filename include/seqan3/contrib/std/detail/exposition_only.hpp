// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::stl::detail implementation helper that are used in multiple files.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_DETAIL_EXPOSITION_ONLY
#define SEQAN_STD_DETAIL_EXPOSITION_ONLY

#include <ranges>

namespace seqan::stl::detail
{

template <typename T>
inline constexpr bool is_initializer_list = false;

template <typename T>
inline constexpr bool is_initializer_list<std::initializer_list<T>> = true;

template <class range_t>
concept simple_view = std::ranges::view<range_t> && std::ranges::range<range_t const>
                   && std::same_as<std::ranges::iterator_t<range_t>, std::ranges::iterator_t<range_t const>>
                   && std::same_as<std::ranges::sentinel_t<range_t>, std::ranges::sentinel_t<range_t const>>;

template <bool is_const, typename t>
using maybe_const = std::conditional_t<is_const, t const, t>;

} // namespace seqan::stl::detail

#endif // SEQAN_STD_DETAIL_EXPOSITION_ONLY
