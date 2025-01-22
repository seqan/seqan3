// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <variant>

#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3
{

/*!\brief A std::variant can be printed by visiting the stream operator for the corresponding type.
 * \tparam variant_ts The types of the variant.
 * \ingroup core_debug_stream
 */
template <typename... variant_ts>
struct std_variant_printer<std::variant<variant_ts...>>
{
    /*!\brief Prints the variant by visiting the stream operator for the corresponding type.
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \param[in,out] stream The output stream.
     * \param[in] arg The variant argument to print.
     *
     * Note that in case the variant is valueless(_by_exception), nothing is printed.
     */
    template <typename stream_t, typename arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        if (!arg.valueless_by_exception())
            std::visit(
                [&stream](auto && elem)
                {
                    stream << elem;
                },
                arg);
        else
            stream << "<VALUELESS_VARIANT>";
    }
};
} // namespace seqan3
