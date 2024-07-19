// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Includes the aligned_sequence and the related insert_gap and
 *        erase_gap functions to enable stl container support.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/range/sequence.hpp>
#include <seqan3/utility/container/concept.hpp>

// ---------------------------------------------------------------------------------------------------------------------
// unaligned_seq transformation trait
// ---------------------------------------------------------------------------------------------------------------------

namespace seqan3::detail
{

//!\brief Helper function to deduce the unaligned sequence type from an aligned sequence container.
template <template <typename...> typename container_type, typename seq_alph_t, typename... rest_t>
    requires container<container_type<gapped<seq_alph_t>, rest_t...>>
constexpr auto
    remove_gap_from_value_type(container_type<gapped<seq_alph_t>, rest_t...>) -> container_type<seq_alph_t, rest_t...>;

//!\overload
template <template <typename...> typename container_type,
          template <typename...>
          typename allocator_type,
          typename seq_alph_t,
          typename... rest_t>
    requires container<container_type<gapped<seq_alph_t>, allocator_type<gapped<seq_alph_t>>, rest_t...>>
constexpr auto
    remove_gap_from_value_type(container_type<gapped<seq_alph_t>, allocator_type<gapped<seq_alph_t>>, rest_t...>)
        -> container_type<seq_alph_t, allocator_type<seq_alph_t>, rest_t...>;

//!\brief Default transformation trait that shall expose the unaligned sequence type of t when specialised.
template <typename t>
struct unaligned_seq
{};

//!\brief Exposes the unaligned sequence type given an aligned sequence container type.
template <typename t>
    requires (!requires { typename std::remove_reference_t<t>::unaligned_sequence_type; })
          && requires { remove_gap_from_value_type(std::declval<t>()); }
struct unaligned_seq<t>
{
    //!\brief The unaligned sequence type of t
    using type = decltype(remove_gap_from_value_type(std::declval<t>()));
};

// customisation point for our gap decorators.
//!\brief Exposes the unaligned sequence type if *t* exposes the type member `unaligned_sequence_type`.
template <typename t>
    requires requires { typename std::remove_reference_t<t>::unaligned_sequence_type; }
struct unaligned_seq<t>
{
    using type = typename std::remove_reference_t<t>::unaligned_sequence_type; //!< The unaligned sequence type of t
};

//!\brief Helper type that delegates to seqan3::detail::unaligned_seq::type.
template <typename t>
using unaligned_seq_t = typename unaligned_seq<t>::type;

} // namespace seqan3::detail

// ---------------------------------------------------------------------------------------------------------------------
// aligned_sequence
// ---------------------------------------------------------------------------------------------------------------------

namespace seqan3
{

/*!\interface seqan3::aligned_sequence <>
 * \extends   seqan3::sequence
 * \brief     The generic concept for an aligned sequence.
 * \ingroup   alignment_aligned_sequence
 *
 * This concept describes the requirements a sequence must fulfil in order to represent a sequence in an alignment.
 *
 * The type needs to model seqan3::sequence and the reference type must be comparable to seqan3::gap.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and type traits.
 * Types that model this concept are shown as "implementing this interface".
 */
/*!\name Requirements for seqan3::aligned_sequence
 * \brief You can expect these functions on all types that model seqan3::aligned_sequence.
 * \relates seqan3::aligned_sequence
 * \{
 */
//!\}
//!\cond
template <typename t>
concept aligned_sequence = sequence<t> && std::equality_comparable_with<std::ranges::range_reference_t<t>, gap>;
//!\endcond

/*!\interface seqan3::writable_aligned_sequence <>
 * \extends   seqan3::aligned_sequence
 * \extends   std::ranges::forward_range
 * \brief     The generic concept for an aligned sequence that is writable.
 * \ingroup   alignment_aligned_sequence
 *
 * This concept describes the requirements a sequence must fulfill in order to be used inside of the alignment algorithm
 * to store the final alignment.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and type traits.
 * Types that model this concept are shown as "implementing this interface".
 */
/*!\name Requirements for seqan3::writable_aligned_sequence
 * \brief You can expect these functions on all types that model seqan3::writable_aligned_sequence.
 * \relates seqan3::writable_aligned_sequence
 * \{
 */
/*!\fn      inline std::ranges::iterator_t<aligned_seq_t> insert_gap(aligned_seq_t & aligned_seq,
 *          typename aligned_seq_t::const_iterator pos_it)
 * \brief   Insert a seqan3::gap into an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::writable_aligned_sequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert a gap.
 * \returns       An iterator pointing to the inserted gap.
 *
 * \details
 * \note      This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 *
 * \attention This is a concept requirement, not an actual function (however types
 *            modeling this concept will provide an implementation).
 */
/*!\fn      inline std::ranges::iterator_t<aligned_seq_t> insert_gap(aligned_seq_t & aligned_seq,
 *          typename aligned_seq_t::const_iterator pos_it, typename aligned_seq_t::size_type size)
 * \brief   Insert multiple seqan3::gap into an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::writable_aligned_sequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert a gaps.
 * \param[in]     size            The number of gap symbols to insert (will result in a gap of length `size`).
 * \returns       An iterator pointing to the first inserted gap or `pos_it` if `size == 0`.
 *
 * \details
 * \note      This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 *
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline std::ranges::iterator_t<aligned_seq_t> erase_gap(aligned_seq_t & aligned_seq,
 *          typename aligned_seq_t::const_iterator pos_it)
 * \brief   Erase a seqan3::gap from an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::writable_aligned_sequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to erase a gap.
 * \returns       An iterator following the removed element. If the iterator `pos_it` refers to the last element, the
 *                std::ranges::end() iterator is returned.
 *
 * \throws seqan3::gap_erase_failure if there is no seqan3::gap at \p pos_it.
 *
 * \details
 * \note      This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 *
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline std::ranges::iterator_t<aligned_seq_t> erase_gap(aligned_seq_t & aligned_seq,
 *          typename aligned_seq_t::const_iterator first, typename aligned_seq_t::const_iterator last)
 * \brief   Erase multiple seqan3::gap from an aligned sequence.
 *
 * \tparam        aligned_seq_t   Type of the range to modify; must model seqan3::writable_aligned_sequence.
 * \param[in,out] aligned_seq     The aligned sequence to modify.
 * \param[in]     first           The iterator pointing to the position where to start erasing gaps.
 * \param[in]     last            The iterator pointing to the position where to stop erasing gaps.
 * \returns       An iterator following the last removed element. If the iterator `last` refers to the last element, the
 *                std::ranges::end() iterator is returned.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) no seqan3::gap.
 *
 * \details
 * \note      This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 *
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      void assign_unaligned(aligned_seq_t & aligned_seq, unaligned_sequence_type && unaligned_seq)
 * \brief   Assign an ungapped sequence to a gapped sequence.
 *
 * \tparam        aligned_seq_t     Type of the container to reassign; must model seqan3::writable_aligned_sequence.
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
//!\cond
template <typename t>
concept writable_aligned_sequence = aligned_sequence<t> && std::ranges::forward_range<t> && requires {
    typename detail::unaligned_seq_t<t>;
} && requires (t v, detail::unaligned_seq_t<t> unaligned) {
    // global functions for generic usability
    { insert_gap(v, std::ranges::begin(v)) } -> std::same_as<std::ranges::iterator_t<t>>;
    { insert_gap(v, std::ranges::begin(v), 2) } -> std::same_as<std::ranges::iterator_t<t>>;
    { erase_gap(v, std::ranges::begin(v)) } -> std::same_as<std::ranges::iterator_t<t>>;
    { erase_gap(v, std::ranges::begin(v), std::ranges::end(v)) } -> std::same_as<std::ranges::iterator_t<t>>;
    { assign_unaligned(v, unaligned) } -> std::same_as<void>;
};
//!\endcond

// ---------------------------------------------------------------------------------------------------------------------
// Aligned sequence interface for containers over the seqan3::gapped alphabet
// ---------------------------------------------------------------------------------------------------------------------

/*!\name Aligned sequence interface for containers over the seqan3::gapped alphabet
 * \brief Enables containers to model seqan3::writable_aligned_sequence if they model seqan3::sequence_container
 *        and have a value type of the seqan3::gapped alphabet.
 * \{
 */
/*!\brief An implementation of seqan3::writable_aligned_sequence::insert_gap for sequence containers.
 * \ingroup alignment_aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::sequence_container;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert a gap.
 * \returns       An iterator pointing to the inserted gap.
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, value)` of
 * the container.
 *
 * \note This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 */
template <sequence_container aligned_seq_t>
    requires detail::is_gapped_alphabet<std::iter_value_t<aligned_seq_t>>
inline typename aligned_seq_t::iterator insert_gap(aligned_seq_t & aligned_seq,
                                                   typename aligned_seq_t::const_iterator pos_it)
{
    return aligned_seq.insert(pos_it, std::iter_value_t<aligned_seq_t>{gap{}});
}

/*!\brief An implementation of seqan3::writable_aligned_sequence::insert_gap for sequence containers.
 * \ingroup alignment_aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::sequence_container;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to insert gaps.
 * \param[in]     size            The number of gap symbols to insert (will result in a gap of length `size`).
 * \returns       An iterator pointing to the first inserted gap or `pos_it` if `size == 0`.
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, `size`, `value`)`
 * of the container.
 *
 * \note This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 */
template <sequence_container aligned_seq_t>
    requires detail::is_gapped_alphabet<std::iter_value_t<aligned_seq_t>>
inline typename aligned_seq_t::iterator insert_gap(aligned_seq_t & aligned_seq,
                                                   typename aligned_seq_t::const_iterator pos_it,
                                                   typename aligned_seq_t::size_type size)
{
    return aligned_seq.insert(pos_it, size, std::iter_value_t<aligned_seq_t>{gap{}});
}

/*!\brief An implementation of seqan3::writable_aligned_sequence::erase_gap for sequence containers.
 * \ingroup alignment_aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::sequence_container;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     pos_it          The iterator pointing to the position where to erase a gap.
 * \returns       An iterator following the removed element. If the iterator `pos_it` refers to the last element, the
 *                std::ranges::end() iterator is returned.
 *
 * \throws seqan3::gap_erase_failure if there is no seqan3::gap at \p pos_it.
 *
 * \details
 *
 * This function delegates to the member function `erase(iterator)` of the
 * container. Before delegating, the function checks if the position pointed to
 * by \p pos_it is an actual seqan3::gap and throws an exception if not.
 *
 * \note This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 */
template <sequence_container aligned_seq_t>
    requires detail::is_gapped_alphabet<std::iter_value_t<aligned_seq_t>>
inline typename aligned_seq_t::iterator erase_gap(aligned_seq_t & aligned_seq,
                                                  typename aligned_seq_t::const_iterator pos_it)
{
    if (*pos_it != gap{}) // [[unlikely]]
        throw gap_erase_failure("The position to be erased does not contain a gap.");

    return aligned_seq.erase(pos_it);
}

/*!\brief An implementation of seqan3::writable_aligned_sequence::erase_gap for sequence containers.
 * \ingroup alignment_aligned_sequence
 * \tparam        aligned_seq_t   Type of the container to modify; must model seqan3::sequence_container;
 *                                The value type must be a seqan3::gapped alphabet.
 * \param[in,out] aligned_seq     The aligned container to modify.
 * \param[in]     first           The iterator pointing to the position where to start erasing gaps.
 * \param[in]     last            The iterator pointing to the position where to stop erasing gaps.
 * \returns       An iterator following the last removed element. If the iterator `last` refers to the last element, the
 *                std::ranges::end() iterator is returned.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) no seqan3::gap.
 *
 * \details
 *
 * This function delegates to the member function `erase(iterator, iterator)` of
 * the container. Before delegating, the function checks if the range
 * [\p first, \p last) contains only seqan3::gap symbols.
 *
 * \note This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 */
template <sequence_container aligned_seq_t>
    requires detail::is_gapped_alphabet<std::iter_value_t<aligned_seq_t>>
inline typename aligned_seq_t::iterator erase_gap(aligned_seq_t & aligned_seq,
                                                  typename aligned_seq_t::const_iterator first,
                                                  typename aligned_seq_t::const_iterator last)
{
    for (auto it = first; it != last; ++it)
        if (*it != gap{}) // [[unlikely]]
            throw gap_erase_failure("The range to be erased contains at least one non-gap character.");

    return aligned_seq.erase(first, last);
}

/*!\brief An implementation of seqan3::writable_aligned_sequence::assign_unaligned_sequence for sequence containers.
 * \ingroup alignment_aligned_sequence
 * \tparam        aligned_seq_t     Type of the container to reassign; must model seqan3::sequence_container;
 *                                  the value type must be a seqan3::gapped alphabet.
 * \tparam        unaligned_seq_t   Type of the container to assign from; must model std::ranges::forward_range.
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
template <sequence_container aligned_seq_t, std::ranges::forward_range unaligned_sequence_type>
    requires detail::is_gapped_alphabet<std::iter_value_t<aligned_seq_t>>
          && weakly_assignable_from<std::ranges::range_reference_t<aligned_seq_t>,
                                    std::ranges::range_reference_t<unaligned_sequence_type>>
inline void assign_unaligned(aligned_seq_t & aligned_seq, unaligned_sequence_type && unaligned_seq)
{
    using std::swap;
    aligned_seq_t tmp;
    tmp.resize(std::ranges::distance(unaligned_seq));
#if SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif // SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
    std::ranges::copy(unaligned_seq, std::ranges::begin(tmp));
#if SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic pop
#endif // SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
    swap(aligned_seq, tmp);
}
//!\}

/*!\name Aligned sequence interface for ranges that have the corresponding member functions.
 * \brief Enables ranges to model seqan3::aligned_sequence if they have the member functions
 *        insert_gap() and erase_gap().
 * \{
 */
/*!\brief An implementation of seqan3::writable_aligned_sequence::insert_gap for ranges with the corresponding
 *        member function insert_gap(it, size).
 * \ingroup alignment_aligned_sequence
 * \tparam range_type    Type of the range to modify; must have an insert_gap(it, size) member function.
 * \param[in,out] rng    The range to modify.
 * \param[in]     pos_it The iterator pointing to the position where to start inserting gaps.
 * \param[in]     size   The number of gaps to insert as an optional argument, default is 1.
 * \returns       An iterator pointing to the first inserted gap or `pos_it` if `size == 0`.
 *
 * \details
 *
 * This function delegates to the member function `insert(iterator, size)` of the range.
 *
 * \note This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 */
template <typename range_type>
    requires requires (range_type v) {
        v.insert_gap(std::ranges::iterator_t<range_type>{});
        v.insert_gap(std::ranges::iterator_t<range_type>{}, typename range_type::size_type{});
    }
std::ranges::iterator_t<range_type> insert_gap(range_type & rng,
                                               std::ranges::iterator_t<range_type> const pos_it,
                                               typename range_type::size_type const size = 1)
{
    return rng.insert_gap(pos_it, size);
}

/*!\brief An implementation of seqan3::writable_aligned_sequence::erase_gap for ranges with the corresponding
 *        member function erase_gap(it).
 * \ingroup alignment_aligned_sequence
 * \tparam range_type   Type of the range to modify; must have an erase_gap(it) member function.
 * \param[in,out] rng   The range to modify.
 * \param[in] pos_it    The iterator pointing to the position where to erase one gap.
 * \returns An iterator following the removed element. If the iterator `pos_it` refers to the last element, the
 *          std::ranges::end() iterator is returned.
 *
 * \details
 *
 * This function delegates to the member function `erase(it)` of
 * the range.
 *
 * \note This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 */
template <typename range_type>
    requires requires (range_type v) { v.erase_gap(std::ranges::iterator_t<range_type>{}); }
std::ranges::iterator_t<range_type> erase_gap(range_type & rng, std::ranges::iterator_t<range_type> const pos_it)
{
    return rng.erase_gap(pos_it);
}

/*!\brief An implementation of seqan3::writable_aligned_sequence::erase_gap for ranges with the corresponding
 *        member function erase_gap(first, last).
 * \ingroup alignment_aligned_sequence
 * \tparam range_type   Type of the range to modify; must have an erase_gap(first, last) member function.
 * \param[in,out] rng   The range to modify.
 * \param[in] first     The iterator pointing to the position where to start erasing gaps.
 * \param[in] last      The iterator pointing to the position where to stop erasing gaps.
 * \returns An iterator following the last removed element. If the iterator `last` refers to the last element, the
 *          std::ranges::end() iterator is returned.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) is no seqan3::gap.
 *
 * \details
 *
 * This function delegates to the member function `erase(first, last)` of the range.
 *
 * \note This may cause reallocations and thus invalidates all iterators and references. Use the returned iterator.
 */
template <typename range_type>
    requires requires (range_type v) {
        v.erase_gap(std::ranges::iterator_t<range_type>{}, std::ranges::iterator_t<range_type>{});
    }
std::ranges::iterator_t<range_type> erase_gap(range_type & rng,
                                              std::ranges::iterator_t<range_type> const first,
                                              std::ranges::iterator_t<range_type> const last)
{
    return rng.erase_gap(first, last);
}
//!\}
} // namespace seqan3

namespace seqan3::detail
{

/*!\brief True, if each type models seqan3::aligned_sequence; false otherwise.
 * \tparam elems The pack of types to be tested.
 */
template <typename... elems>
inline constexpr bool all_model_aligned_seq = (aligned_sequence<elems> && ...);

/*!\brief True, if each type in the seqan3::type_list models seqan3::aligned_sequence; false otherwise.
 * \tparam elems The pack of types to be tested.
 */
template <typename... elems>
inline constexpr bool all_model_aligned_seq<type_list<elems...>> = all_model_aligned_seq<elems...>;
} // namespace seqan3::detail
