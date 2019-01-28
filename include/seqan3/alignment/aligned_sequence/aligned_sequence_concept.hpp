// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Includes the aligned_sequence_concept and the related insert_gap and
 *        erase_gap functions to enable stl container support.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iomanip>
#include <tuple>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/slice.hpp>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

// -----------------------------------------------------------------------------
// aligned_sequence_concept
// -----------------------------------------------------------------------------

/*!\interface seqan3::aligned_sequence_concept <>
 * \extends   std::ranges::ForwardRange
 * \brief     The generic concept for an aligned sequence.
 * \ingroup   aligned_sequence
 *
 * This concept describes the requirements a sequence must fulfil
 * in order to be part of the seqan3::alignment object.
 *
 * The following extended type requirements for a type `T` must hold true:
 *
 *   * seqan3::reference_t<T> must model seqan3::alphabet_concept.
 *   * seqan3::reference_t<T> must be assignable from seqan3::gap.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that model this concept are shown as "implementing this interface".
 */
/*!\name Requirements for seqan3::aligned_sequence_concept
 * \brief You can expect these functions on all types that model seqan3::aligned_sequence_concept.
 * \relates seqan3::aligned_sequence_concept
 * \{
 */
/*!\fn      inline seq_type::iterator insert_gap(seq_type & seq, typename seq_type::const_iterator pos_it)
 * \brief   Insert a seqan3::gap into an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to insert a gap.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline seq_type::iterator insert_gap(seq_type & seq, typename seq_type::const_iterator pos_it, typename seq_type::size_type size)
 * \brief   Insert multiple seqan3::gap into an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to insert a gaps.
 * \param[in]     size     The number of gap symbols to insert (will result in a gap of length `size`).
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline seq_type::iterator erase_gap(seq_type & seq, typename seq_type::const_iterator pos_it)
 * \brief   Erase a seqan3::gap from an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to erase a gap.
 *
 * \throws seqan3::gap_erase_failure if there is no seqan3::gap at \p pos_it.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline seq_type::iterator erase_gap(seq_type & seq, typename seq_type::const_iterator first, typename seq_type::const_iterator last)
 * \brief   Erase multiple seqan3::gap from an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     first    The iterator pointing to the position where to start erasing gaps.
 * \param[in]     last     The iterator pointing to the position where to stop erasing gaps.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) no seqan3::gap.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
//!\}
//!\cond
template <typename t>
SEQAN3_CONCEPT aligned_sequence_concept =
    std::ranges::ForwardRange<t> &&
    alphabet_concept<value_type_t<t>> &&
    weakly_assignable_concept<reference_t<t>, gap const &> &&
    requires (t v)
    {
        { insert_gap(v, v.begin()) } -> typename t::iterator; // global functions for generic usability
        { insert_gap(v, v.begin(), 2) } -> typename t::iterator;
        { erase_gap(v, v.begin()) } -> typename t::iterator;
        { erase_gap(v, v.begin(), v.end()) } -> typename t::iterator;
    };
//!\endcond

// -----------------------------------------------------------------------------
// Functions that make sequence containers model aligned_sequence_concept
// -----------------------------------------------------------------------------

/*!\name Aligned sequence interface for containers
 * \brief Enables containers to model seqan3::aligned_sequence_concept if they
 * model seqan3::sequence_container_concept and have a reference type assignable
 * from seqan3::gap (e.g. std::vector<seqan3::gapped<seqan3::dna4>>).
 * \{
 */
/*!\brief An implementation of seqan3::aligned_sequence_concept::insert_gap for sequence containers.
 * \tparam        seq_type Type of the container to modify; must model
 *                         seqan3::sequence_container_concept; the reference type
 *                         (seqan3::reference_t<seq_type>) must be assignable from
 *                         seqan3::gap.
 * \param[in,out] seq      The container to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to insert a gap.
 *
 * \relates seqan3::aligned_sequence_concept
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, value)` of
 * the container.
 */
template <sequence_container_concept seq_type>
//!\cond
    requires weakly_assignable_concept<reference_t<seq_type>, gap const &>
//!\endcond
inline typename seq_type::iterator insert_gap(seq_type & seq, typename seq_type::const_iterator pos_it)
{
    return seq.insert(pos_it, value_type_t<seq_type>{gap{}});
}

/*!\brief An implementation of seqan3::aligned_sequence_concept::insert_gap for sequence containers.
 * \tparam        seq_type Type of the container to modify; must model
 *                         seqan3::sequence_container_concept; the reference type
 *                         (seqan3::reference_t<seq_type>) must be assignable from
 *                         seqan3::gap.
 * \param[in,out] seq      The container to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to insert gaps.
 * \param[in]     size     The number of gap symbols to insert (will result in a gap of length `size`).
 *
 * \relates seqan3::aligned_sequence_concept
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, `size`, `value`)`
 * of the container.
 */
template <sequence_container_concept seq_type>
//!\cond
    requires weakly_assignable_concept<reference_t<seq_type>, gap const &>
//!\endcond
inline typename seq_type::iterator insert_gap(seq_type & seq,
                                              typename seq_type::const_iterator pos_it,
                                              typename seq_type::size_type size)
{
    return seq.insert(pos_it, size, value_type_t<seq_type>{gap{}});
}

/*!\brief An implementation of seqan3::aligned_sequence_concept::erase_gap for sequence containers.
 * \tparam        seq_type Type of the container to modify; must model
 *                         seqan3::sequence_container_concept; the reference type
 *                         (seqan3::reference_t<seq_type>) must be assignable from
 *                         seqan3::gap.
 * \param[in,out] seq      The container to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to erase a gap.
 *
 * \throws seqan3::gap_erase_failure if there is no seqan3::gap at \p pos_it.
 *
 * \relates seqan3::aligned_sequence_concept
 *
 * \details
 *
 * This function delegates to the member function `erase(iterator)` of the
 * container. Before delegating, the function checks if the position pointed to
 * by \p pos_it is an actual seqan3::gap and throws an exception if not.
 */
template <sequence_container_concept seq_type>
//!\cond
    requires weakly_assignable_concept<reference_t<seq_type>, gap const &>
//!\endcond
inline typename seq_type::iterator erase_gap(seq_type & seq, typename seq_type::const_iterator pos_it)
{
    if (*pos_it != gap{}) // [[unlikely]]
        throw gap_erase_failure("The position to be erased does not contain a gap.");

    return seq.erase(pos_it);
}

/*!\brief An implementation of seqan3::aligned_sequence_concept::erase_gap for sequence containers.
 * \tparam        seq_type Type of the container to modify; must model
 *                         seqan3::sequence_container_concept; the reference type
 *                         (seqan3::reference_t<seq_type>) must be assignable from
 *                         seqan3::gap.
 * \param[in,out] seq      The container to modify.
 * \param[in]     first    The iterator pointing to the position where to start erasing gaps.
 * \param[in]     last     The iterator pointing to the position where to stop erasing gaps.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) no seqan3::gap.
 *
 * \relates seqan3::aligned_sequence_concept
 *
 * \details
 *
 * This function delegates to the member function `erase(iterator, iterator)` of
 * the container. Before delegating, the function checks if the range
 * [\p first, \p last) contains only seqan3::gap symbols.
 */
template <sequence_container_concept seq_type>
//!\cond
    requires weakly_assignable_concept<reference_t<seq_type>, gap const &>
//!\endcond
inline typename seq_type::iterator erase_gap(seq_type & seq,
                                             typename seq_type::const_iterator first,
                                             typename seq_type::const_iterator last)
{
    for (auto it = first; it != last; ++it)
        if (*it != gap{}) // [[unlikely]]
            throw gap_erase_failure("The range to be erased contains at least one non-gap character.");

    return seq.erase(first, last);
}
//!\}

namespace detail
{

/*!\brief               Create the formatted alignment output and add it to the provided debug_stream.
 * \ingroup             aligned_sequence
 * \tparam alignment_t  The type of the alignment, must satisfy tuple_like_concept.
 * \tparam idx          An index sequence.
 * \param[in] stream    The output stream that receives the formatted alignment.
 * \param[in] align     The alignment that shall be streamed.
 */
template<tuple_like_concept alignment_t, size_t ...idx>
void stream_alignment(debug_stream_type & stream, alignment_t const & align, std::index_sequence<idx...> const & /**/)
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
        ranges::for_each(get<0>(align) | ranges::view::slice(begin_pos, end_pos) | view::to_char,
                         [&stream] (char ch) { stream << ch; });

        auto stream_f = [&] (auto const & previous_seq, auto const & aligned_seq)
        {
            // write alignment bars
            stream << '\n' << std::setw(8) << "";
            ranges::for_each(ranges::zip_view(previous_seq, aligned_seq) | ranges::view::slice(begin_pos, end_pos),
                             [&stream] (auto && ch) { stream << (get<0>(ch) == get<1>(ch) ? '|' : ' '); });

            // write next sequence
            stream << '\n' << std::setw(8) << "";
            ranges::for_each(aligned_seq | ranges::view::slice(begin_pos, end_pos) | view::to_char,
                             [&stream] (char ch) { stream << ch; });
        };
        (stream_f(get<idx>(align), get<idx + 1>(align)), ...);
        stream << '\n';
    }
}

/*!\brief True, if each type satisfies aligned_sequence_concept; false otherwise.
 * \tparam elems The pack of types to be tested.
 */
template <typename ...elems>
inline bool constexpr all_satisfy_aligned_seq = false;

/*!\brief True, if each type satisfies aligned_sequence_concept; false otherwise.
 * \tparam elems The pack of types to be tested.
 */
template <typename ...elems>
inline bool constexpr all_satisfy_aligned_seq<type_list<elems...>> = (aligned_sequence_concept<elems> && ...);

} // namespace detail

/*!\brief Streaming operator for alignments, which are represented as tuples of aligned sequences.
 * \tparam tuple_t  The alignment type, must satisfy tuple_like_concept and its size must be at least 2.
 * \param stream    The target stream for the formatted output.
 * \param alignment The alignment that shall be formatted. All sequences must be equally long.
 * \return          The given stream to which the alignment representation is appended.
 */
template <tuple_like_concept tuple_t>
//!\cond
    requires detail::all_satisfy_aligned_seq<detail::tuple_type_list_t<tuple_t>>
//!\endcond
inline debug_stream_type & operator<<(debug_stream_type & stream, tuple_t const & alignment)
{
    static_assert(std::tuple_size_v<tuple_t> >= 2, "An alignment requires at least two sequences.");
    detail::stream_alignment(stream, alignment, std::make_index_sequence<std::tuple_size_v<tuple_t> - 1> {});
    return stream;
}

} // namespace seqan
