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
#include <cstdio>

#include "../../alphabet/nucleotide/dna4_container.hpp"

namespace seqan3
{

template <typename t>
concept bool aligned_sequence_concept = true;

//! Alignment class.
/*!
    An alignment is a tuple of at least two aligned sequences.
*/
template <typename ...sequence_ts>
    requires (aligned_sequence_concept<sequence_ts> && ...)
class alignment : public std::tuple<sequence_ts...>
{
public:
    //! The alignment depth is the number of sequences contained.
    std::size_t const depth = sizeof...(sequence_ts);

    //! Constructor that allows to pass sequences directly.
    constexpr alignment(sequence_ts && ...sequences)
        : std::tuple<sequence_ts...>(std::forward<sequence_ts>(sequences)...)
    {
        static_assert(sizeof...(sequence_ts) > 1, "An alignment requires at least two sequences.");
    };
};

//! Alignment class template deduction as for tuple.
template <typename ...sequence_ts>
alignment(sequence_ts...) -> alignment<sequence_ts...>;

namespace detail
{
//! Create the formatted alignment output and add it to a stream.
template <typename stream_t, typename tuple_t, std::size_t ...idx>
void stream_alignment(stream_t & stream, tuple_t const & tuple, std::index_sequence<idx...> const & /**/)
{
    std::size_t const alignment_length = std::get<0>(tuple).size();
    char buf[9];

    // split alignment into blocks of length 50 and loop over parts
    for (std::size_t used_length = 0; used_length < alignment_length; used_length += 50)
    {
        // write header
        std::snprintf(buf, 9, "%7lu ", used_length);
        stream << std::endl << buf;
        for (std::size_t col = 1; col <= 50 && col + used_length <= alignment_length; ++col)
        {
            if (col % 10 == 0)
                stream << ':';
            else if (col % 5 == 0)
                stream << '.';
            else
                stream << ' ';
        }

        // write sequences
        const char * indent = "        ";
        stream << std::endl << indent << std::get<0>(tuple).substr(used_length, 50);
        auto stream_f = [&]
            (auto const & previous_sequence, auto const & aligned_sequence)
        {
            stream << std::endl << indent;
            auto seq1 = previous_sequence.begin() + used_length;
            auto seq2 = aligned_sequence.begin() + used_length;
            for (auto it1 = seq1, it2 = seq2;
                 it1 < previous_sequence.end() && it1 < seq1 + 50 &&
                 it2 < aligned_sequence.end()  && it2 < seq2 + 50;
                 ++it1, ++it2)
            {
                stream << (*it1 == *it2 ? '|' : ' ');
            }
            stream << std::endl << indent << aligned_sequence.substr(used_length, 50);
        };
        (stream_f(std::get<idx>(tuple), std::get<idx + 1>(tuple)), ...);
        stream << std::endl;
    }
}

} // namespace detail

//! Formatted output stream of an alignment.
template <typename stream_type, typename ...sequence_ts>
stream_type & operator<<(stream_type & outstream, alignment<sequence_ts...> const & align)
{
    static_assert(sizeof...(sequence_ts) >= 2, "An alignment requires at least two sequences.");
    detail::stream_alignment(outstream, align, std::make_index_sequence<sizeof...(sequence_ts) - 1>{});
    return outstream;
}

} // namespace seqan3
