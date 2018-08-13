// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

/*!\file
 * \brief Provides the seqan3::alignment class.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

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

#include <seqan3/range/view/to_char.hpp>

namespace seqan3
{

// TODO(joergi-w) Remove this when aligned_sequence_concept is ready.
template <typename t>
concept bool aligned_sequence_concept = true;

/*!
 * \brief An alignment is a tuple of at least two aligned sequences.
 * \tparam sequence_types Types of the aligned sequences. Each must satisfy the seqan3::aligned_sequence_concept.
 */
template <typename ...sequence_types>
    requires (aligned_sequence_concept<sequence_types> && ...)
class alignment : public std::tuple<sequence_types...>
{
public:
    //! The alignment depth is the number of sequences contained.
    std::size_t const depth = sizeof...(sequence_types);

    //! Constructor that allows to pass sequences directly.
    constexpr alignment(sequence_types && ...sequences)
        : std::tuple<sequence_types...>(std::forward<sequence_types>(sequences)...)
    {
        static_assert(sizeof...(sequence_types) > 1, "An alignment requires at least two sequences.");
    };

    // Constructors, destructor and assignment
    alignment(alignment const &) = default;
    alignment & operator=(alignment const &) = default;
    alignment(alignment &&) = default;
    alignment & operator=(alignment &&) = default;
};

/*!
 * Alignment class template deduction.
 * \tparam sequence_types Types of the aligned sequences. Each must satisfy the seqan3::aligned_sequence_concept.
 * \sa http://en.cppreference.com/w/cpp/language/class_template_deduction
 */
template <typename ...sequence_types>
alignment(sequence_types...) -> alignment<sequence_types...>;

//! Type for column_view.
template <typename ...sequence_types>
using column_view_type = ranges::zip_view<ranges::iterator_range<typename sequence_types::const_iterator,
                                                                 typename sequence_types::const_iterator>...>;

/*!
 * \brief Column-wise view of a sequence alignment.
 * \tparam sequence_types Types of the aligned sequences. Each must satisfy the seqan3::aligned_sequence_concept.
 * \param align The underlying alignment for the view.
 * \return A view containing the alignment columns for the given alignment.
 *
 * ```cpp
 *     alignment align("AUUGN"_rna5, "AGUGN"_rna5);
 *     for (std::tuple<rna5, rna5> const & col : column_view(align))
 *     {
 *         // col -> {A, A}, {U, G}, ...
 *     }
 * ```
 */
template <typename ...sequence_types>
    requires (aligned_sequence_concept<sequence_types> && ...)
column_view_type<sequence_types...> column_view(alignment<sequence_types...> const & align)
{
    return std::apply([] (auto && ...args) { return ranges::zip_view(ranges::view::all(args)...); },
                      align);
}

namespace detail
{
/*!
 * Create the formatted alignment output and add it to a stream.
 * \tparam stream_t Type of the output stream.
 * \tparam alignment_t Type of the alignment.
 * \tparam idx Index sequence.
 * \param stream The output stream that receives the formatted alignment.
 * \param align The alignment that shall be streamed.
 */
template <typename stream_t, typename alignment_t, std::size_t ...idx>
void stream_alignment(stream_t & stream, alignment_t const & align, std::index_sequence<idx...> const & /**/)
{
    std::size_t const alignment_length = std::get<0>(align).size();

    // split alignment into blocks of length 50 and loop over parts
    for (std::size_t used_length = 0; used_length < alignment_length; used_length += 50)
    {
        // write header
        if (used_length != 0)
            stream << std::endl;

        stream << std::setw(7) << used_length << ' ';
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
                stream << (seqan3::to_char(*it1) == seqan3::to_char(*it2) ? '|' : ' ');
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
 * \tparam stream_type    Type of the output stream.
 * \tparam sequence_types Types of the aligned sequences.
 * \param  outstream      The output stream that receives the formatted alignment.
 * \param  alignment      The alignment that shall be streamed.
 * \return                The output stream that contains the formatted alignment.
 */
template <typename stream_type, typename ...sequence_types>
stream_type & operator<<(stream_type & outstream, alignment<sequence_types...> const & alignment)
{
    static_assert(sizeof...(sequence_types) >= 2, "An alignment requires at least two sequences.");
    detail::stream_alignment(outstream, alignment, std::make_index_sequence<sizeof...(sequence_types) - 1> {});
    return outstream;
}

} // namespace seqan3

namespace std
{

/*!
 * tuple_size overload for alignment class (to enable functions like std::apply).
 * \tparam sequence_types Types of the aligned sequences.
 */
template <typename ...sequence_types>
struct tuple_size<seqan3::alignment<sequence_types...>> : tuple_size<tuple<sequence_types...>> {};

} // namespace std
