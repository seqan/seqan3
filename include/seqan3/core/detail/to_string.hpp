// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Auxiliary for pretty printing of exception messages.
 */

#pragma once

#include <seqan3/core/debug_stream.hpp>

namespace seqan3::detail
{

/*!\brief Streams all parameters via the seqan3::debug_stream and returns a concatenated string.
 *
 * \tparam    value_type Must be streamable (stream << value).
 * \param[in] values     Variable number of parameters of any type that implement the stream operator.
 * \returns A concatenated string of all values (no separator in between is added).
 */
template <typename ...value_type>
std::string to_string(value_type && ...values)
{
    std::stringstream stream;
    debug_stream_type dstream{stream};
    (dstream << ... << std::forward<value_type>(values));
    return stream.str();
}

} // namespace seqan3::detail
