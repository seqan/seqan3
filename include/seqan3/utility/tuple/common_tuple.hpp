// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
using common_tuple = seqan::std::tuple<t...>;

} // namespace seqan3
