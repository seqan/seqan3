// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Includes the AlignedSequence and the related insert_gap and
 *        erase_gap functions to enable stl container support.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iomanip>
#include <tuple>

#include <range/v3/algorithm/for_each.hpp>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/type_traits/all.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/slice.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>

// ---------------------------------------------------------------------------------------------------------------------
// unaligned_seq transformation trait
// ---------------------------------------------------------------------------------------------------------------------

namespace seqan3::detail
{

//!\brief Helper function to deduce the unaligned sequence type from an aligned sequence container.
template <template <typename...> typename container_type, typename seq_alph_t, typename ...rest_t>
//!\cond
    requires Container<container_type<gapped<seq_alph_t>, rest_t...>>
//!\endcond
constexpr auto remove_gap_from_value_type(container_type<gapped<seq_alph_t>, rest_t...>)
    -> container_type<seq_alph_t, rest_t...>;

//!\overload
template <template <typename...> typename container_type,
          template <typename...> typename allocator_type,
          typename seq_alph_t, typename ...rest_t>
//!\cond
    requires Container<container_type<gapped<seq_alph_t>, allocator_type<gapped<seq_alph_t>>, rest_t...>>
//!\endcond
constexpr auto remove_gap_from_value_type(container_type<gapped<seq_alph_t>, allocator_type<gapped<seq_alph_t>>, rest_t...>)
    -> container_type<seq_alph_t, allocator_type<seq_alph_t>, rest_t...>;

//!\brief Default transformation trait that shall expose the unaligned sequence type of t when specialised.
template <typename t>
struct unaligned_seq
{};

//!\brief Exposes the unaligned sequence type given an aligned sequence container type.
template <typename t>
//!\cond
    requires !requires { typename std::remove_reference_t<t>::unaligned_seq_type; } &&
              requires { remove_gap_from_value_type(std::declval<t>()); }
//!\endcond
struct unaligned_seq<t>
{
    //!\brief The unaligned sequence type of t
    using type = decltype(remove_gap_from_value_type(std::declval<t>()));
};

// customisation point for our gap decorators.
//!\brief Exposes the unaligned sequence type if *t* exposes the type member `unaligned_seq_type`.
template <typename t>
//!\cond
    requires requires { typename std::remove_reference_t<t>::unaligned_seq_type; }
//!\endcond
struct unaligned_seq<t>
{
    using type = typename std::remove_reference_t<t>::unaligned_seq_type; //!< The unaligned sequence type of t
};

//!\brief Helper type that delegates to seqan3::detail::unaligned_seq::type.
template <typename t>
using unaligned_seq_t = typename unaligned_seq<t>::type;

} // namespace seqan3::detail

// ---------------------------------------------------------------------------------------------------------------------
// AlignedSequence
// ---------------------------------------------------------------------------------------------------------------------

namespace seqan3
{

/*!\interface seqan3::AlignedSequence <>
 * \extends   std::ranges::ForwardRange
 * \brief     The generic concept for an aligned sequence.
 * \ingroup   aligned_sequence
 *
 * This concept describes the requirements a sequence must fulfil in order to be used inside of the alignment algorithm
 * to store the final alignment.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and type traits.
 * Types that model this concept are shown as "implementing this interface".
 */
/*!\name Requirements for seqan3::AlignedSequence
 * \brief You can expect these functions on all types that model seqan3::AlignedSequence.
 * \relates seqan3::AlignedSequence
 * \{
 */
/*!\fn      inline aligned_seq_t::iterator insert_gap(aligned_seq_t & aligned_seq, typename aligned_seq_t::const_iterator pos_it)
 * \brief   Insert a seqan3::gap into an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::AlignedSequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert a gap.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline aligned_seq_t::iterator insert_gap(aligned_seq_t & aligned_seq, typename aligned_seq_t::const_iterator
 *          pos_it, typename aligned_seq_t::size_type size)
 * \brief   Insert multiple seqan3::gap into an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::AlignedSequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert a gaps.
 * \param[in]     size            The number of gap symbols to insert (will result in a gap of length `size`).
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline aligned_seq_t::iterator erase_gap(aligned_seq_t & aligned_seq, typename aligned_seq_t::const_iterator
 *          pos_it)
 * \brief   Erase a seqan3::gap from an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::AlignedSequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to erase a gap.
 *
 * \throws seqan3::gap_erase_failure if there is no seqan3::gap at \p pos_it.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline aligned_seq_t::iterator erase_gap(aligned_seq_t & aligned_seq, typename aligned_seq_t::const_iterator
 *          first, typename aligned_seq_t::const_iterator last)
 * \brief   Erase multiple seqan3::gap from an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::AlignedSequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     first           The iterator pointing to the position where to start erasing gaps.
 * \param[in]     last            The iterator pointing to the position where to stop erasing gaps.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) no seqan3::gap.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      void assign_unaligned(aligned_seq_t & aligned_seq, unaligned_sequence_type && unaligned_seq)
 * \brief   Assign an ungapped sequence to a gapped sequence.
 *
 * \tparam        aligned_seq_t     Type of the container to reassign; must model seqan3::AlignedSequence.
 * \tparam        unaligned_seq_t   Type of the container to assign from; must correspond to the aligned type without
 *                                  gap information (see details.)
 * \param[in,out] aligned_seq       The gapped sequence container to assign to.
 * \param[in,out] unaligned_seq     The unaligned sequence container to assign from.
 *
 * \details
 *
 * An aligned sequence has to be assignable from its unaligned counter part. For example a
 * std::vector<seqan3::gapped<seqan3::dna4>> as well as a seqan3::gap_decorator<std::vector<seqan3::dna4>>
 * can be assigned from s std::vector<seqan3::dna4> via seqan3::assign_unaligned.
 *
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
//!\}
//!\cond
template <typename t>
SEQAN3_CONCEPT AlignedSequence =
    std::ranges::ForwardRange<t> &&
    Alphabet<reference_t<t>> &&
    WeaklyAssignable<reference_t<t>, gap const &> &&
    requires { typename detail::unaligned_seq_t<t>; } &&
    requires (t v, detail::unaligned_seq_t<t> unaligned)
    {
        { insert_gap(v, v.begin()) } -> typename t::iterator; // global functions for generic usability
        { insert_gap(v, v.begin(), 2) } -> typename t::iterator;
        { erase_gap(v, v.begin()) } -> typename t::iterator;
        { erase_gap(v, v.begin(), v.end()) } -> typename t::iterator;
        { assign_unaligned(v, unaligned) } -> void;
    };
//!\endcond

// ---------------------------------------------------------------------------------------------------------------------
// Aligned sequence interface for containers over the seqan3::gapped alphabet
// ---------------------------------------------------------------------------------------------------------------------

/*!\name Aligned sequence interface for containers over the seqan3::gapped alphabet
 * \brief Enables containers to model seqan3::AlignedSequence if they model seqan3::SequenceContainer
 *        and have a value type of the seqan3::gapped alphabet.
 * \{
 */
/*!\brief An implementation of seqan3::AlignedSequence::insert_gap for sequence containers.
 * \ingroup aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::SequenceContainer;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert a gap.
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, value)` of
 * the container.
 */
template <SequenceContainer aligned_seq_t>
//!\cond
    requires detail::is_gapped_alphabet<value_type_t<aligned_seq_t>>
//!\endcond
inline typename aligned_seq_t::iterator insert_gap(aligned_seq_t & aligned_seq,
                                                   typename aligned_seq_t::const_iterator pos_it)
{
    return aligned_seq.insert(pos_it, value_type_t<aligned_seq_t>{gap{}});
}

/*!\brief An implementation of seqan3::AlignedSequence::insert_gap for sequence containers.
 * \ingroup aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::SequenceContainer;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert gaps.
 * \param[in]     size            The number of gap symbols to insert (will result in a gap of length `size`).
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, `size`, `value`)`
 * of the container.
 */
template <SequenceContainer aligned_seq_t>
//!\cond
    requires detail::is_gapped_alphabet<value_type_t<aligned_seq_t>>
//!\endcond
inline typename aligned_seq_t::iterator insert_gap(aligned_seq_t & aligned_seq,
                                                   typename aligned_seq_t::const_iterator pos_it,
                                                   typename aligned_seq_t::size_type size)
{
    return aligned_seq.insert(pos_it, size, value_type_t<aligned_seq_t>{gap{}});
}

/*!\brief An implementation of seqan3::AlignedSequence::erase_gap for sequence containers.
 * \ingroup aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::SequenceContainer;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to erase a gap.
 *
 * \throws seqan3::gap_erase_failure if there is no seqan3::gap at \p pos_it.
 *
 * \details
 *
 * This function delegates to the member function `erase(iterator)` of the
 * container. Before delegating, the function checks if the position pointed to
 * by \p pos_it is an actual seqan3::gap and throws an exception if not.
 */
template <SequenceContainer aligned_seq_t>
//!\cond
    requires detail::is_gapped_alphabet<value_type_t<aligned_seq_t>>
//!\endcond
inline typename aligned_seq_t::iterator erase_gap(aligned_seq_t & aligned_seq,
                                                  typename aligned_seq_t::const_iterator pos_it)
{
    if (*pos_it != gap{}) // [[unlikely]]
        throw gap_erase_failure("The position to be erased does not contain a gap.");

    return aligned_seq.erase(pos_it);
}

/*!\brief An implementation of seqan3::AlignedSequence::erase_gap for sequence containers.
 * \ingroup aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::SequenceContainer;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     first           The iterator pointing to the position where to start erasing gaps.
 * \param[in]     last            The iterator pointing to the position where to stop erasing gaps.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) no seqan3::gap.
 *
 * \details
 *
 * This function delegates to the member function `erase(iterator, iterator)` of
 * the container. Before delegating, the function checks if the range
 * [\p first, \p last) contains only seqan3::gap symbols.
 */
template <SequenceContainer aligned_seq_t>
//!\cond
    requires detail::is_gapped_alphabet<value_type_t<aligned_seq_t>>
//!\endcond
inline typename aligned_seq_t::iterator erase_gap(aligned_seq_t & aligned_seq,
                                                  typename aligned_seq_t::const_iterator first,
                                                  typename aligned_seq_t::const_iterator last)
{
    for (auto it = first; it != last; ++it)
        if (*it != gap{}) // [[unlikely]]
            throw gap_erase_failure("The range to be erased contains at least one non-gap character.");

    return aligned_seq.erase(first, last);
}

/*!\brief An implementation of seqan3::AlignedSequence::assign_unaligned_sequence for sequence containers.
 * \ingroup aligned_sequence
 * \tparam        aligned_seq_t     Type of the container to reassign; must model seqan3::SequenceContainer;
 *                                  the value type must be a seqan3::gapped alphabet.
 * \tparam        unaligned_seq_t   Type of the container to assign from; must model std::ranges::ForwardRange.
 * \param[in,out] aligned_seq       The gapped sequence container to assign to.
 * \param[in,out] unaligned_seq     The unaligned sequence container to assign from.
 *
 * \details
 *
 * This function clears the content of the `gapped` container and reassigns the content of the `unaligned` container by
 * using std::copy.
 *
 * ### Performance
 *
 * Linear in the size of unaligned_seq.
 *
 * ### Exceptions
 *
 * Strong exception guarantee.
 *
 */
template <SequenceContainer aligned_seq_t, std::ranges::ForwardRange unaligned_sequence_type>
//!\cond
    requires detail::is_gapped_alphabet<value_type_t<aligned_seq_t>> &&
             WeaklyAssignable<reference_t<aligned_seq_t>, reference_t<unaligned_sequence_type>>
//!\endcond
inline void assign_unaligned(aligned_seq_t & aligned_seq, unaligned_sequence_type && unaligned_seq)
{
    using std::swap;
    aligned_seq_t tmp;
    tmp.resize(std::ranges::distance(unaligned_seq));
    std::copy(std::ranges::begin(unaligned_seq), std::ranges::end(unaligned_seq), std::ranges::begin(tmp));
    swap(aligned_seq, tmp);
}
//!\}

/*!\name Aligned sequence interface for ranges that have the corresponding member functions.
 * \brief Enables ranges to model seqan3::AlignedSequence if they have the member functions
 *        insert_gap() and erase_gap().
 * \{
 */
/*!\brief An implementation of seqan3::AlignedSequence::insert_gap for ranges with the corresponding
 *        member function insert_gap(it, size).
 * \ingroup aligned_sequence
 * \tparam range_type   Type of the range to modify; must have an insert_gap(it, size) member function.
 * \param[in,out] rng   The range to modify.
 * \param[in]     it    The iterator pointing to the position where to start inserting gaps.
 * \param[in]     size  The number of gaps to insert as an optional argument, default is 1.
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, size)` of the range.
 */
template <typename range_type>
//!\cond
    requires requires (range_type v)
        {
            v.insert_gap(typename range_type::iterator{});
            v.insert_gap(typename range_type::iterator{}, typename range_type::size_type{});
        }
//!\endcond
typename range_type::iterator insert_gap(range_type & rng,
                                         typename range_type::iterator const it,
                                         typename range_type::size_type const size = 1)
{
    return rng.insert_gap(it, size);
}

/*!\brief An implementation of seqan3::AlignedSequence::erase_gap for ranges with the corresponding
 *        member function erase_gap(it).
 * \ingroup aligned_sequence
 * \tparam range_type   Type of the range to modify; must have an erase_gap(it) member function.
 * \param[in,out] rng   The range to modify.
 * \param[in] it        The iterator pointing to the position where to erase one gap.
 *
 * \details
 *
 * This function delegates to the member function `erase(it)` of
 * the range.
 */
template <typename range_type>
//!\cond
    requires requires (range_type v) { v.erase_gap(typename range_type::iterator{}); }
//!\endcond
typename range_type::iterator erase_gap(range_type & rng,
                                        typename range_type::iterator const it)
{
    return rng.erase_gap(it);
}

/*!\brief An implementation of seqan3::AlignedSequence::erase_gap for ranges with the corresponding
 *        member function erase_gap(first, last).
 * \ingroup aligned_sequence
 * \tparam range_type   Type of the range to modify; must have an erase_gap(first, last) member function.
 * \param[in,out] rng   The range to modify.
 * \param[in] first     The iterator pointing to the position where to start erasing gaps.
 * \param[in] last      The iterator pointing to the position where to stop erasing gaps.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) is no seqan3::gap.
 *
 * \details
 *
 * This function delegates to the member function `erase(first, last)` of the range.
 */
template <typename range_type>
//!\cond
    requires requires (range_type v) { v.erase_gap(typename range_type::iterator{}, typename range_type::iterator{}); }
//!\endcond
typename range_type::iterator erase_gap(range_type & rng,
                                        typename range_type::iterator const first,
                                        typename range_type::iterator const last)
{
    return rng.erase_gap(first, last);
}
//!\}

namespace detail
{

/*!\brief               Create the formatted alignment output and add it to the provided debug_stream.
 * \ingroup             aligned_sequence
 * \tparam alignment_t  The type of the alignment, must satisfy TupleLike.
 * \tparam idx          An index sequence.
 * \param[in] stream    The output stream that receives the formatted alignment.
 * \param[in] align     The alignment that shall be streamed.
 */
template<TupleLike alignment_t, size_t ...idx>
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
        ranges::for_each(get<0>(align) | view::slice(begin_pos, end_pos) | view::to_char,
                         [&stream] (char ch) { stream << ch; });

        auto stream_f = [&] (auto const & previous_seq, auto const & aligned_seq)
        {
            // write alignment bars
            stream << '\n' << std::setw(8) << "";
            ranges::for_each(ranges::zip_view(previous_seq, aligned_seq) | view::slice(begin_pos, end_pos),
                             [&stream] (auto && ch) { stream << (get<0>(ch) == get<1>(ch) ? '|' : ' '); });

            // write next sequence
            stream << '\n' << std::setw(8) << "";
            ranges::for_each(aligned_seq | view::slice(begin_pos, end_pos) | view::to_char,
                             [&stream] (char ch) { stream << ch; });
        };
        (stream_f(get<idx>(align), get<idx + 1>(align)), ...);
        stream << '\n';
    }
}

/*!\brief True, if each type satisfies AlignedSequence; false otherwise.
 * \tparam elems The pack of types to be tested.
 */
template <typename ...elems>
inline bool constexpr all_satisfy_aligned_seq = false;

/*!\brief True, if each type satisfies AlignedSequence; false otherwise.
 * \tparam elems The pack of types to be tested.
 */
template <typename ...elems>
inline bool constexpr all_satisfy_aligned_seq<type_list<elems...>> = (AlignedSequence<elems> && ...);

} // namespace detail

/*!\brief Streaming operator for alignments, which are represented as tuples of aligned sequences.
 * \ingroup         aligned_sequence
 * \tparam tuple_t  The alignment type, must satisfy TupleLike and its size must be at least 2.
 * \param stream    The target stream for the formatted output.
 * \param alignment The alignment that shall be formatted. All sequences must be equally long.
 * \return          The given stream to which the alignment representation is appended.
 */
template <TupleLike tuple_t>
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
