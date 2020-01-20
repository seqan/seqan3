// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides pretty printing for chars.
 */

#pragma once

#include <string>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// make_printable
// ----------------------------------------------------------------------------

/*!\brief Returns a printable value for the given character `c`.
 * \param[in] c The character to be represented as printable string.
 * \return    a std::string containing a printable version of the given character `c`.
 * \ingroup char_operations
 *
 * \details
 *
 * Some characters, e.g. control commands, cannot be printed. This function converts them to a std::string
 * containing the visual representation of this character. For all control commands the value `'CTRL'` is returned.
 *
 * ### Exception
 *
 * Strong exception guarantee is given.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
inline std::string make_printable(char const c)
{
    switch (c)
    {
        case '\0':                   return "'\\0'";
        case '\t':                   return "'\\t'";
        case '\n':                   return "'\\n'";
        case '\v':                   return "'\\v'";
        case '\f':                   return "'\\f'";
        case '\r':                   return "'\\r'";
        case static_cast<char>(127): return "'DEL'";
        default:
        {
            if ((c >= static_cast<char>(1) && c <= static_cast<char>(8)) ||
                (c >= static_cast<char>(14) && c <= static_cast<char>(31)))
                return "'CTRL'";
            else
                return {'\'', c, '\''};
        }
    }
}

} // namespace seqan3::detail
