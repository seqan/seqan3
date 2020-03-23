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

#include <iostream>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/detail/debug_stream_optional.hpp>
#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/core/detail/debug_stream_tuple.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/detail/debug_stream_variant.hpp>

// forward declare
//!\cond
namespace std
{
extern ostream cerr;
} // namespace std
//!\endcond

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::debug_stream
// ------------------------------------------------------------------

//!\brief A global instance of seqan3::debug_stream_type.
//!\ingroup stream
inline debug_stream_type debug_stream{std::cerr};

} // namespace seqan3
