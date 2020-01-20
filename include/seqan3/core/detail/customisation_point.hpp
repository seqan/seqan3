// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Helper utilities for defining customisation point objects.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ============================================================================
// priority_tag
// ============================================================================

//!\brief A tag that allows controlled overload resolution via implicit base conversion rules.
template <size_t I>
struct priority_tag
//!\cond
// Doxygen fail
: priority_tag<I-1>
//!\endcond
{};

//!\brief Recursion anchor for seqan3::detail::priority_tag.
template <>
struct priority_tag<0>
{};

} // seqan3::detail

// ============================================================================
// SEQAN3_CPO_IMPL
// ============================================================================

//!\brief A macro that helps defining the overload set of a customisation point.
#define SEQAN3_CPO_IMPL(PRIO, TERM)                                                                                  \
/*!\brief A customisation point overload.*/                                                                          \
template <typename t, typename ...arg_ts>                                                                            \
static constexpr decltype(auto) impl(seqan3::detail::priority_tag<PRIO>,                                             \
                                     [[maybe_unused]] t && v,                                                        \
                                     [[maybe_unused]] arg_ts && ... args)                                            \
    noexcept(noexcept(TERM))                                                                                         \
    requires requires (seqan3::detail::priority_tag<PRIO> const &/*<- need for doxygen*/, t && v, arg_ts && ... args)\
    { { TERM }; }                                                                                                    \
{                                                                                                                    \
    return TERM;                                                                                                     \
}
