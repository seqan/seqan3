// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#pragma once

/*! \file
 *  \ingroup io
 *  \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

 #include <fstream>
 #include <sstream>

namespace seqan3
{
// basic meta_functions:
template <typename t>
concept bool stream_concept = requires(t stream)
{
    typename t::char_type;
    typename t::traits_type;
    typename t::int_type;
    typename t::pos_type;
    typename t::off_type;

    // Requires to access underlying stream buffer.
    { stream.rdbuf() } -> std::basic_streambuf<typename t::char_type, typename t::traits_type>*;
};

template <typename t>
concept bool input_stream_concept = stream_concept<t> &&
    requires(t stream, typename t::char_type val)
{
    // formatted input.
    { operator>>(stream, val) } -> auto&;
};

template <typename t>
concept bool output_stream_concept = stream_concept<t> &&
    requires(t stream)
{
    // formatted output.
    { operator<<(stream, typename t::char_type{}) } -> auto&;
};

}  // namespace seqan3

//#ifndef NDEBUG

namespace seqan3::detail
{
#include <type_traits>
#include <istream>
#include <ostream>
#include <vector>

// check input_stream_concept
static_assert(input_stream_concept<std::istream>);
static_assert(!input_stream_concept<std::ostream>);
static_assert(stream_concept<std::istream>);
static_assert(input_stream_concept<std::iostream>);
static_assert(stream_concept<std::iostream>);
static_assert(!input_stream_concept<std::vector<char>>);

// check output_stream_concept
static_assert(!output_stream_concept<std::istream>);
static_assert(output_stream_concept<std::ostream>);
static_assert(stream_concept<std::ostream>);
static_assert(output_stream_concept<std::iostream>);
static_assert(stream_concept<std::iostream>);
static_assert(!output_stream_concept<std::vector<char>>);
}

//#endif // NDEBUG
