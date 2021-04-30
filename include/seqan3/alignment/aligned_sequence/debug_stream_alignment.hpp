// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief The seqan3::debug_stream_type overload in order to print alignments.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iomanip>
#include <tuple>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/detail/debug_stream_tuple.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief               Create the formatted alignment output and add it to the provided debug_stream.
 * \ingroup             aligned_sequence
 * \tparam alignment_t  The type of the alignment; must model seqan3::tuple_like.
 * \tparam idx          An index sequence.
 * \param[in] stream    The output stream that receives the formatted alignment.
 * \param[in] align     The alignment that shall be streamed.
 */
template <typename char_t, tuple_like alignment_t, size_t ...idx>
void stream_alignment(debug_stream_type<char_t> & stream, alignment_t const & align, std::index_sequence<idx...> const & /**/)
{
    using std::get;
    size_t const alignment_size = get<0>(align).size();

    // split alignment into blocks of length 50 and loop over parts
    for (size_t begin_pos = 0; begin_pos < alignment_size; begin_pos += 50)
    {
        size_t const end_pos = std::min(begin_pos + 50, alignment_size);

        // write header line
        if (begin_pos != 0)
            stream << '\n';

        stream << std::setw(7) << begin_pos << ' ';
        for (size_t pos = begin_pos + 1; pos <= end_pos; ++pos)
        {
            if (pos % 10 == 0)
                stream << ':';
            else if (pos % 5 == 0)
                stream << '.';
            else
                stream << ' ';
        }

        // write first sequence
        stream << '\n' << std::setw(8) << "";
        std::ranges::for_each(get<0>(align) | views::slice(begin_pos, end_pos) | views::to_char,
                              [&stream] (char ch) { stream << ch; });

        auto stream_f = [&] (auto const & previous_seq, auto const & aligned_seq)
        {
            // write alignment bars
            stream << '\n' << std::setw(8) << "";
            std::ranges::for_each(views::zip(previous_seq, aligned_seq) | views::slice(begin_pos, end_pos),
                                  [&stream] (auto && ch) { stream << (get<0>(ch) == get<1>(ch) ? '|' : ' '); });

            // write next sequence
            stream << '\n' << std::setw(8) << "";
            std::ranges::for_each(aligned_seq | views::slice(begin_pos, end_pos) | views::to_char,
                                  [&stream] (char ch) { stream << ch; });
        };
        (stream_f(get<idx>(align), get<idx + 1>(align)), ...);
        stream << '\n';
    }
}
} // namespace seqan3::detail

namespace seqan3
{
/*!\brief Stream operator for alignments, which are represented as tuples of aligned sequences.
 * \ingroup aligned_sequence
 *
 * \tparam alignment_t The alignment type, must satisfy tuple_like and its size must be at least 2.
 *
 * \param[in,out] stream The target stream for the formatted output.
 * \param[in] alignment The alignment that shall be formatted. All sequences must be equally long.
 *
 * \return The given stream to which the alignment representation is appended.
 */
template <typename char_t, typename alignment_t>
//!\cond
    requires (detail::debug_streamable_tuple<alignment_t> &&
              detail::all_model_aligned_seq<detail::tuple_type_list_t<std::remove_cvref_t<alignment_t>>>)
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & stream, alignment_t && alignment)
{
    constexpr size_t sequence_count = std::tuple_size_v<std::remove_cvref_t<alignment_t>>;

    static_assert(sequence_count >= 2, "An alignment requires at least two sequences.");

    detail::stream_alignment(stream, alignment, std::make_index_sequence<sequence_count - 1>{});
    return stream;
}

} // namespace seqan
