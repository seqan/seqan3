// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides literals for unit conversions.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::test
{

inline namespace literals
{

//!\brief A string literal that interprets the given number as MiB and converts to bytes.
static constexpr size_t operator""_MiB(unsigned long long int const number) noexcept
{
    return number << 20;
}

//!\brief A string literal that interprets the given number as GiB and converts to bytes.
static constexpr size_t operator""_GiB(unsigned long long int const number) noexcept
{
    return number << 30;
}

} // namespace literals

} // namespace seqan3::test
