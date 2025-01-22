// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <concepts>
#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
} // namespace seqan3

namespace seqan3::detail
{

//!\brief Helper function to print elements of a tuple separately.
template <typename char_t, typename tuple_t, std::size_t... I>
void print_tuple(debug_stream_type<char_t> & s, tuple_t && t, std::index_sequence<I...> const &)
{
    using std::get;
    s << '(';
    ((s << (I == 0 ? "" : ",") << get<I>(t)), ...);
    s << ')';
}

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Printer for formatted output of tuple like objects.
 * \tparam tuple_t Type of the tuple to be printed; must model seqan3::tuple_like.
 * \ingroup core_debug_stream
 */
template <tuple_like tuple_t>
struct tuple_printer<tuple_t>
{
    /*!\brief Prints a tuple to a formatted output stream.
     *
     * Takes a stream and a tuple object as arguments and prints the tuple elements to the stream.
     * The elements are printed using the `detail::print_tuple` function.
     *
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \param[in,out] stream The stream to print to.
     * \param[in] arg The tuple to print.
     */
    template <typename stream_t, typename arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        detail::print_tuple(stream, std::forward<arg_t>(arg), std::make_index_sequence<std::tuple_size_v<tuple_t>>{});
    }
};

//!\}

} // namespace seqan3
