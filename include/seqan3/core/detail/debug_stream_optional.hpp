// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <optional>

#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
/*!\brief Make std::nullopt_t printable.
 * \tparam    optional_type This is std::nullopt_t.
 * \param[in] s             The seqan3::debug_stream.
 * \param[in] arg           This is std::nullopt.
 * \relates seqan3::debug_stream_type
 */
template <typename optional_type, typename char_t>
//!\cond
   requires std::same_as<remove_cvref_t<optional_type>, std::nullopt_t>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, optional_type && SEQAN3_DOXYGEN_ONLY(arg))
{
   s << "<VALUELESS_OPTIONAL>";
   return s;
}

/*!\brief A std::optional can be printed by printing its value or nothing if valueless.
 * \tparam    optional_type The type of the optional.
 * \param[in] s             The seqan3::debug_stream.
 * \param[in] arg           The std::optional.
 * \relates seqan3::debug_stream_type
 */
template <typename optional_type, typename char_t>
//!\cond
    requires detail::is_type_specialisation_of_v<remove_cvref_t<optional_type>, std::optional>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, optional_type && arg)
{
    if (arg.has_value())
        s << *arg;
    else
        s << "<VALUELESS_OPTIONAL>";
    return s;
}

//!\}

} // namespace seqan3
