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

#include <iomanip>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/iterator_range.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/view/slice.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet.hpp>
#include <seqan3/range/view/to_char.hpp>

namespace seqan3
{

template <typename t>
concept bool aligned_sequence_concept = true;

//! Alignment class.
/*!
 * An alignment is a tuple of at least two aligned sequences.
 * \tparam sequence_ts sequence types of the alignment
 */
template <typename ...sequence_ts>
    requires (aligned_sequence_concept<sequence_ts> && ...)
class alignment : public std::tuple<sequence_ts...>
{
public:
    //! The alignment depth is the number of sequences contained.
    static std::size_t const depth = sizeof...(sequence_ts);

    //! Constructor that allows to pass sequences directly.
    constexpr alignment(sequence_ts && ...sequences)
        : std::tuple<sequence_ts...>(std::forward<sequence_ts>(sequences)...)
    {
        static_assert(sizeof...(sequence_ts) > 1, "An alignment requires at least two sequences.");
    };
};

/*!
 * Alignment class template deduction.
 * \tparam sequence_ts sequence types of the alignment
 * \sa http://en.cppreference.com/w/cpp/language/class_template_deduction
 */
template <typename ...sequence_ts>
alignment(sequence_ts...) -> alignment<sequence_ts...>;

//! Type for column_iterator.
template <typename ...sequence_ts>
using column_iterator_type = ranges::v3::zip_view<ranges::v3::iterator_range<typename sequence_ts::const_iterator,
                                                                             typename sequence_ts::const_iterator>...>;
/*!
 * Column-wise iteration over a sequence alignment.
 * \tparam sequence_ts sequence types of the alignment
 * \param align the underlying alignment for the iterator
 * \return a column iterator for the given alignment
 */
template <typename ...sequence_ts>
    requires (aligned_sequence_concept<sequence_ts> && ...)
column_iterator_type<sequence_ts...> column_iterator(alignment<sequence_ts...> const & align)
{
    return std::apply([] (auto && ...args) { return ranges::v3::zip_view(ranges::v3::view::all(args)...); },
                      align);
}

namespace detail
{
/*!
 * Create the formatted alignment output and add it to a stream.
 * \tparam stream_t output stream type
 * \tparam alignment_t alignment type
 * \tparam idx index sequence
 * \param stream output stream that receives the formatted alignment
 * \param align the alignment that shall be streamed
 */
template <typename stream_t, typename alignment_t, std::size_t ...idx>
void stream_alignment(stream_t & stream, alignment_t const & align, std::index_sequence<idx...> const & /**/)
{
    std::size_t const alignment_length = std::get<0>(align).size();

    // split alignment into blocks of length 50 and loop over parts
    for (std::size_t used_length = 0; used_length < alignment_length; used_length += 50)
    {
        // write header
        stream << std::endl << std::setw(7) << used_length << ' ';
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
        stream << std::endl << indent;
        std::size_t const col_end = std::min(used_length + 50, alignment_length);
        ranges::for_each(std::get<0>(align) | ranges::view::slice(used_length, col_end) | view::to_char,
                         [&stream] (char ch) { stream << ch; });

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
            stream << std::endl << indent;
            ranges::for_each(aligned_sequence | ranges::view::slice(used_length, col_end) | view::to_char,
                             [&stream] (char ch) { stream << ch; });
        };
        (stream_f(std::get<idx>(align), std::get<idx + 1>(align)), ...);
        stream << std::endl;
    }
}

} // namespace detail

/*!
 * Formatted output stream of an alignment.
 * \tparam stream_type output stream type
 * \tparam sequence_ts sequence types of the alignment
 * \param outstream output stream that receives the formatted alignment
 * \param align the alignment that shall be streamed
 * \return output stream that contains the formatted alignment
 */
template <typename stream_type, typename ...sequence_ts>
stream_type & operator<<(stream_type & outstream, alignment<sequence_ts...> const & align)
{
    static_assert(sizeof...(sequence_ts) >= 2, "An alignment requires at least two sequences.");
    detail::stream_alignment(outstream, align, std::make_index_sequence<sizeof...(sequence_ts) - 1> {});
    return outstream;
}

} // namespace seqan3

namespace std
{

/*!
 * tuple_size overload for alignment class (to enable functions like std::apply).
 * \tparam sequence_ts sequence types of the alignment
 */
template <typename ...sequence_ts>
struct tuple_size<seqan3::alignment<sequence_ts...>> : tuple_size<tuple<sequence_ts...>> {};

} // namespace std
