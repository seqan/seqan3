// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <iostream>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/debug_stream/all.hpp>

// forward declare
//!\cond
namespace std
{
namespace
{
extern ostream cerr;
}
} // namespace std
//!\endcond

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::debug_stream
// ------------------------------------------------------------------

//!\brief A global instance of seqan3::debug_stream_type.
//!\ingroup core_debug_stream
inline debug_stream_type debug_stream{std::cerr};

} // namespace seqan3
