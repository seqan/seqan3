// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various shortcuts for common std::ranges functions.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

//!\cond
template <typename t>
SEQAN3_DEPRECATED_310 auto begin(t && v)
{
    return std::ranges::begin(std::forward<t>(v));
}

template <typename t>
SEQAN3_DEPRECATED_310 auto end(t && v)
{
    return std::ranges::end(std::forward<t>(v));
}

template <typename t>
SEQAN3_DEPRECATED_310 auto cbegin(t && v)
{
    return std::ranges::cbegin(std::forward<t>(v));
}

template <typename t>
SEQAN3_DEPRECATED_310 auto cend(t && v)
{
    return std::ranges::cend(std::forward<t>(v));
}

template <typename t>
SEQAN3_DEPRECATED_310 auto size(t && v)
{
    return std::ranges::size(std::forward<t>(v));
}

template <typename t>
SEQAN3_DEPRECATED_310 auto empty(t && v)
{
    return std::ranges::empty(std::forward<t>(v));
}
//!\endcond

} // namespace seqan3
