// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
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
/*!\name Formatted output overloads
 * \{
 */

template <typename T>
    requires (!std::is_same_v<T, std::remove_cvref_t<T>>)
struct optional_printer<T> : public optional_printer<std::remove_cvref_t<T>>
{
};

/*!\brief Make std::nullopt_t printable.
 * \tparam    optional_type This is std::nullopt_t.
 * \param[in] s             The seqan3::debug_stream.
 * \param[in] arg           This is std::nullopt.
 * \relates seqan3::debug_stream_type
 */
template <>
struct optional_printer<std::nullopt_t>
{
    constexpr static auto print = [](auto & s, std::nullopt_t SEQAN3_DOXYGEN_ONLY(arg))
    {
        s << "<VALUELESS_OPTIONAL>";
    };
};

/*!\brief A std::optional can be printed by printing its value or nothing if valueless.
 * \tparam    optional_type The type of the optional.
 * \param[in] s             The seqan3::debug_stream.
 * \param[in] arg           The std::optional.
 * \relates seqan3::debug_stream_type
 */
template <typename T>
struct optional_printer<std::optional<T>>
{
    constexpr static auto print = [](auto & s, auto && arg)
    {
        if (arg.has_value())
            s << *arg;
        else
            s << "<VALUELESS_OPTIONAL>";
    };
};

//!\}

} // namespace seqan3
