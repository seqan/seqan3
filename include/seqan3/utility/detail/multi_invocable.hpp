// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides a type that combines multiple invocables.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief A type that can conveniently inherit multiple invocables and acts as a union over them.
 * \ingroup utility
 * \tparam invocable_ts The types to inherit from.
 */
template <typename... invocable_ts>
struct multi_invocable : invocable_ts...
{
    //!\brief Inherit the function call operators.
    using invocable_ts::operator()...;
};

//!\brief Deduction guides for seqan3::detail::multi_invocable.
template <typename... invocable_ts>
multi_invocable(invocable_ts...) -> multi_invocable<invocable_ts...>;

} // namespace seqan3::detail
