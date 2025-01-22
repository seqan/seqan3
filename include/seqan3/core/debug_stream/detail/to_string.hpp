// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Auxiliary for pretty printing of exception messages.
 */

#pragma once

#include <sstream>

#include <seqan3/core/debug_stream/debug_stream_type.hpp>

namespace seqan3::detail
{

/*!\brief Streams all parameters via the seqan3::debug_stream and returns a concatenated string.
 *
 * \tparam    value_type Must be streamable (stream << value).
 * \param[in] values     Variable number of parameters of any type that implement the stream operator.
 * \returns A concatenated string of all values (no separator in between is added).
 */
template <typename... value_type>
std::string to_string(value_type &&... values)
{
    std::stringstream stream;
    debug_stream_type dstream{stream};
    (dstream << ... << std::forward<value_type>(values));
    return stream.str();
}

} // namespace seqan3::detail
