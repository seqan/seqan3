// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides various transformation traits for use on iterators.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
/*!\brief Defines iterator_category member if underlying_iterator_t has a valid std::iterator_traits::iterator_category
 *        overload.
 * \ingroup core
 * \details
 *
 * The C++ paper [P2259R1](https://wg21.link/p2259r1) defines the behaviour of iterator_category for a C++-20 input
 * iterator.
 *
 * https://wg21.link/p2259r1#fixing-counted_iterator states:
 *
 * > Provide member iterator_concept and iterator_category when the wrapped iterator type provides them, to honor its
 * > opt-in/opt-outs.
 *
 * That means, only define iterator_category if the underlying iterator has it.
 *
 * \see https://github.com/seqan/product_backlog/issues/151
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96070
 */
template <typename underlying_iterator_t>
struct maybe_iterator_category
{
#if SEQAN3_DOXYGEN_ONLY(1) 0
    /*!\brief The iterator category tag. (not always present!)
     * \details
     *
     * This member is only defined if and only if std::iterator_traits<underlying_iterator_t>::iterator_category is
     * valid and denotes a type.
     */
    using iterator_category = MAYBE_PRESENT(std::iterator_traits<underlying_iterator_t>::iterator_category);
#endif // SEQAN3_DOXYGEN_ONLY(1)0
};

//!\cond
template <typename t>
concept has_iterator_category = requires () { typename t::iterator_category; };
//!\endcond

//!\cond
template <typename underlying_iterator_t>
    requires has_iterator_category<std::iterator_traits<underlying_iterator_t>>
struct maybe_iterator_category<underlying_iterator_t>
{
    using iterator_category = typename std::iterator_traits<underlying_iterator_t>::iterator_category;
};
//!\endcond

/*!\brief This handles more cases than maybe_iterator_category if you inherit the underling_iterator_t.
 * \ingroup core
 * \details
 *
 * The same as maybe_iterator_category, but this will not define the iterator_category if the underling_iterator_t
 * already defines an iterator_category to avoid duplicated definitions of iterator_category from two different base
 * classes.
 *
 * This class can be safely inherited from.
 */
template <typename underling_iterator_t>
struct maybe_inherited_iterator_category : public maybe_iterator_category<underling_iterator_t>
{};

//!\cond
template <typename underling_iterator_t>
    requires has_iterator_category<underling_iterator_t>
struct maybe_inherited_iterator_category<underling_iterator_t>
{
    // underling_iterator_t::iterator_category is already defined
};
//!\endcond

/*!\brief Exposes the
 * [iterator_concept](https://en.cppreference.com/w/cpp/iterator/iterator_tags) from the modelled concept.
 * \ingroup core
 * \implements seqan3::transformation_trait
 * \tparam it_t The type to operate on.
 */
template <typename it_t>
    requires std::input_or_output_iterator<it_t>
using iterator_concept_tag_t = std::conditional_t<
    std::contiguous_iterator<it_t>,
    std::contiguous_iterator_tag,
    std::conditional_t<std::random_access_iterator<it_t>,
                       std::random_access_iterator_tag,
                       std::conditional_t<std::bidirectional_iterator<it_t>,
                                          std::bidirectional_iterator_tag,
                                          std::conditional_t<std::forward_iterator<it_t>,
                                                             std::forward_iterator_tag,
                                                             std::conditional_t<std::input_iterator<it_t>,
                                                                                std::input_iterator_tag,
                                                                                std::output_iterator_tag>>>>>;

} // namespace seqan3::detail

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// iter_pointer
// ----------------------------------------------------------------------------

/*!\brief This is like std::iter_value_t, but for the pointer type.
 * \ingroup core
 * \implements seqan3::transformation_trait
 * \tparam it_t The type to operate on.
 * \see seqan3::detail::iter_pointer_t
 *
 * \attention
 * C++20 does not provide a `std::iter_pointer_t`, because the new C++20 iterators do not need to provide a pointer
 * type.
 */
template <typename it_t>
struct iter_pointer
{
    //!\brief The pointer type of std::iterator_traits or void.
    using type = void;
};

//!\cond
template <typename it_t>
    requires requires { typename std::iterator_traits<it_t>::pointer; }
struct iter_pointer<it_t>
{
    //!\brief This is defined for every legacy input-iterator.
    //!\sa https://en.cppreference.com/w/cpp/iterator/iterator_traits
    using type = typename std::iterator_traits<it_t>::pointer;
};
//!\endcond

/*!\brief Return the `pointer` type of the input type (transformation_trait shortcut).
 * \ingroup core
 * \tparam it_t The type to operate on.
 * \see seqan3::detail::iter_pointer
 */
template <typename it_t>
using iter_pointer_t = typename iter_pointer<it_t>::type;

} // namespace seqan3::detail
