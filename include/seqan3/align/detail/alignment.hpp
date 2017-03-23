// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================
// Author: Joerg Winkler <j.winkler AT fu-berlin.de>
// ============================================================================

#pragma once

#include <tuple>
#include <type_traits>
#include <iostream>
#include <utility>

namespace seqan3
{

template <typename t>
concept bool aligned_sequence_concept = true;

//! Alignment class.
/*!
    An alignment is a tuple of aligned sequences.
*/
template <typename first_t, typename ...remaining_ts>
    requires aligned_sequence_concept<first_t> &&
             (aligned_sequence_concept<remaining_ts> && ...)
class alignment : public std::tuple<first_t, remaining_ts...>
{
public:
    //! The alignment depth is the number of sequences contained.
    std::size_t const depth = sizeof...(remaining_ts) + 1ul;

    //! Constructor that allows to pass sequences directly.
    constexpr alignment(first_t && v1, remaining_ts && ...v2) : std::tuple<first_t, remaining_ts...>(v1, std::forward<remaining_ts>(v2)...){};
};

//! Alignment class template deduction as for tuple.
template <typename ...ts>
alignment(ts...) -> alignment<ts...>;

namespace detail
{
// Helper for output stream operator.
template <typename stream_t, typename tuple_t, size_t ...idx>
void stream_alignment(stream_t & stream, tuple_t const & tuple, std::index_sequence<idx...> const & /**/)
{
    auto stream_f = [&stream] (auto const & aligned_sequence) { stream << aligned_sequence << '\n'; };
    (stream_f(std::get<idx>(tuple)), ...);
}

} // namespace detail

//! Formatted output stream of an alignment.
template <typename stream_type, typename ...ts>
stream_type & operator<<(stream_type & outstream, alignment<ts...> const & align)
{
    detail::stream_alignment(outstream, align, std::make_index_sequence<sizeof...(ts)>{});
    return outstream;
}

} // namespace seqan3
