// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <optional>

#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3
{

//!\brief Printer for formatted output of std::nullopt_t.
//!\ingroup core_debug_stream
template <>
struct optional_printer<std::nullopt_t>
{
    /*!\brief Prints `std::nullopt_t` to formatted output stream.
     * \tparam stream_t The type of the stream to which the value is streamed.
     * \param[in,out] stream The stream to print to.
     * \param[in] arg The optional to print.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, std::nullopt_t const SEQAN3_DOXYGEN_ONLY(arg)) const
    {
        stream << "<VALUELESS_OPTIONAL>";
    }
};

/*!\brief Printer for formatted output of a std::optional
 * \tparam T The value type of the std::optional.
 * \ingroup core_debug_stream
 */
template <typename T>
struct optional_printer<std::optional<T>>
{
    /*!\brief Print the optional to the stream by printing its value or nothing if valueless.
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \param[in,out] stream The stream to print to.
     * \param[in] arg The optional to print.
     */
    template <typename stream_t, typename arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        if (arg.has_value())
            stream << *arg;
        else
            stream << "<VALUELESS_OPTIONAL>";
    }
};

} // namespace seqan3
