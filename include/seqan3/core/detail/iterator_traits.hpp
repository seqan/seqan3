// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various transformation traits for use on iterators.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/iterator>
#include <type_traits>

#include <seqan3/core/platform.hpp>
#include <seqan3/utility/type_traits/pre.hpp>

namespace seqan3
{

/*!\addtogroup core
 * \{
 */

// ----------------------------------------------------------------------------
// value_type
// ----------------------------------------------------------------------------

#ifdef SEQAN3_DEPRECATED_310
namespace detail
{
/*!\brief Exposes the `value_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::input_iterator.
 * \deprecated This is deprecated use std::iter_value_t.
 */
template <std::input_iterator it_t>
struct value_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_value_t<it_t>;
};
} // namespace seqan3::detail
#endif // SEQAN3_DEPRECATED_310

// see specialisation for ranges in core/range/type_traits.hpp

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

#ifdef SEQAN3_DEPRECATED_310
namespace detail
{
/*!\brief Exposes the `reference` type of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::input_iterator.
 * \deprecated This is deprecated use std::iter_reference_t.
 */
template <std::input_iterator it_t>
struct reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_reference_t<it_t>;
};
} // namespace seqan3::detail
#endif // SEQAN3_DEPRECATED_310

// see specialisation for ranges in core/range/type_traits.hpp

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

#ifdef SEQAN3_DEPRECATED_310
namespace detail
{
/*!\brief Exposes the `rvalue_reference` type of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::input_iterator.
 * \deprecated This is deprecated use std::iter_rvalue_reference_t.
 */
template <std::input_iterator it_t>
struct rvalue_reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_rvalue_reference_t<it_t>;
};
} // namespace seqan3::detail
#endif // SEQAN3_DEPRECATED_310

// see specialisation for ranges in core/range/type_traits.hpp

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

// only defined for ranges

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

#ifdef SEQAN3_DEPRECATED_310
namespace detail
{
/*!\brief Exposes the `difference_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::weakly_incrementable.
 * \deprecated This is deprecated use std::iter_difference_t.
 */
template <std::weakly_incrementable it_t>
struct difference_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_difference_t<it_t>;
};
} // namespace seqan3::detail
#endif // SEQAN3_DEPRECATED_310

// see specialisation for ranges in core/range/type_traits.hpp

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

#ifdef SEQAN3_DEPRECATED_310
namespace detail
{
/*!\brief Exposes the `size_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::weakly_incrementable.
 * \deprecated This is deprecated. There is no alternative! Unlike std::ranges::range_size_t, the Standard has no
 *            std::iter_size_t. We decided that it does not make sense to define it on the difference type of the
 *            iterator.
 */
template <std::weakly_incrementable it_t>
struct size_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::make_unsigned_t<std::iter_difference_t<it_t>>;
};
} // namespace seqan3::detail
#endif // SEQAN3_DEPRECATED_310

// see specialisation for ranges in core/range/type_traits.hpp
//!\}

} // namespace seqan3

namespace seqan3::detail
{
/*!\brief Defines iterator_category member if underlying_iterator_t has a valid std::iterator_traits::iterator_category
 *        overload.
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
#if SEQAN3_DOXYGEN_ONLY(1)0
    /*!\brief The iterator category tag. (not always present!)
     * \details
     *
     * This member is only defined if and only if std::iterator_­traits<underlying_iterator_t>::​iterator_­category is
     * valid and denotes a type.
     */
    using iterator_category = MAYBE_PRESENT(std::iterator_­traits<underlying_iterator_t>::​iterator_­category);
#endif // SEQAN3_DOXYGEN_ONLY(1)0
};

//!\cond
template <typename t>
SEQAN3_CONCEPT has_iterator_category = requires ()
{
    typename t::iterator_category;
};
//!\endcond

#if SEQAN3_WORKAROUND_GCC_96070
/*!\brief This is a workaround for gcc 10.x, x < 4. That version of the stdlib always expects an iterator_category
 *        to be defined. There are some view combinations that do not work with this "fix".
 */
template <typename underlying_iterator_t>
    requires (!has_iterator_category<std::iterator_traits<underlying_iterator_t>>)
struct maybe_iterator_category<underlying_iterator_t>
{
    using iterator_category = void;
};
#endif // SEQAN3_WORKAROUND_GCC_96070

//!\cond
template <typename underlying_iterator_t>
    requires has_iterator_category<std::iterator_traits<underlying_iterator_t>>
struct maybe_iterator_category<underlying_iterator_t>
{
    using iterator_category = typename std::iterator_traits<underlying_iterator_t>::iterator_category;
};
//!\endcond

/*!\brief This handles more cases than maybe_iterator_category if you inherit the underling_iterator_t.
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
 * \implements seqan3::transformation_trait
 * \tparam it_t The type to operate on.
 */
template <typename it_t>
//!\cond
    requires std::input_or_output_iterator<it_t>
//!\endcond
using iterator_concept_tag_t =
    std::conditional_t<
        std::contiguous_iterator<it_t>,
        std::contiguous_iterator_tag,
        std::conditional_t<
            std::random_access_iterator<it_t>,
            std::random_access_iterator_tag,
            std::conditional_t<
                std::bidirectional_iterator<it_t>,
                std::bidirectional_iterator_tag,
                std::conditional_t<
                    std::forward_iterator<it_t>,
                    std::forward_iterator_tag,
                    std::conditional_t<
                        std::input_iterator<it_t>,
                        std::input_iterator_tag,
                        std::output_iterator_tag>>>>>;

} // namespace seqan3::detail

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// iter_pointer
// ----------------------------------------------------------------------------

/*!\brief This is like std::iter_value_t, but for the pointer type.
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
 * \tparam it_t The type to operate on.
 * \see seqan3::detail::iter_pointer
 */
template <typename it_t>
using iter_pointer_t = typename iter_pointer<it_t>::type;

} // namespace seqan3::detail
