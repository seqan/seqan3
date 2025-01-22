// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::common_tuple.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/contrib/std/tuple.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/utility/tuple/common_pair.hpp> // Included to not break API.

namespace seqan3
{

/*!\brief A [std::tuple](https://en.cppreference.com/w/cpp/utility/tuple) implementation that incorporates most changes
 *        from C++23's standard library.
 * \ingroup utility_tuple
 */
template <typename... t>
using common_tuple = seqan::stl::tuple<t...>;

} // namespace seqan3
