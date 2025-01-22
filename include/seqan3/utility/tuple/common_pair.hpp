// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::common_pair.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/contrib/std/pair.hpp>
#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief A [std::pair](https://en.cppreference.com/w/cpp/utility/pair) implementation that incorporates most changes
 *        from C++23's standard library.
 * \ingroup utility_tuple
 */
template <typename t1, typename t2>
using common_pair = seqan::stl::pair<t1, t2>;

} // namespace seqan3
