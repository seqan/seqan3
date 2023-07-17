// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
using common_pair = seqan::std::pair<t1, t2>;

} // namespace seqan3
