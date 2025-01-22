// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides an hacky overload for the googletest PrintTo function.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <ios>
#include <iterator>
#include <ostream>

#include <seqan3/core/debug_stream/debug_stream_type.hpp>

namespace seqan3
{

//!\cond DEV
//!\brief Overload for the googletest PrintTo function that always delegates to our debug_stream.
template <typename t>
    requires (!std::input_or_output_iterator<t> && !std::same_as<t, std::default_sentinel_t>)
// tricks the compiler to consider this as more specialized than googletests generic PrintTo
void PrintTo(t const & v, std::ostream * out)
{
    debug_stream_type my_stream{*out};
    my_stream << v;
}
//!\endcond

} // namespace seqan3

namespace std
{

//!\brief Overload for the googletest PrintTo function that always delegates to our debug_stream.
using ::seqan3::PrintTo;

} // namespace std

namespace seqan3::detail
{

//!\brief Overload for the googletest PrintTo function that always delegates to our debug_stream.
using ::seqan3::PrintTo;

} // namespace seqan3::detail
