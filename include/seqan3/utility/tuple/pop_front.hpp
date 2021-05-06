// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::tuple_pop_front.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/utility/tuple/split.hpp>

namespace seqan3
{
/*!\brief Removes the first element of a tuple.
 * \ingroup utility
 *
 * \param[in] t  The original tuple.
 *
 * \returns A new tuple without the first element of `t`.
 *
 * \details
 *
 * Note, that the tuple must contain at least one element and must support empty tuple types, i.e. std::pair cannot
 * be used.
 *
 * ### Complexity
 *
 * Linear in the number of elements.
 *
 * ### Thread safety
 *
 * Concurrent invocations of this functions are thread safe.
 */
template <tuple_like tuple_t>
constexpr auto tuple_pop_front(tuple_t && t)
{
    static_assert(std::tuple_size_v<std::remove_cvref_t<tuple_t>> > 0);

    return std::get<1>(tuple_split<1>(std::forward<tuple_t>(t)));
}

} // namespace seqan3
