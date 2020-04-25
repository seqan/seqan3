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

#include <iterator>
#include <type_traits>

#include <seqan3/core/platform.hpp>
#include <seqan3/core/type_traits/pre.hpp>
#include <seqan3/std/iterator>

namespace seqan3
{

/*!\addtogroup type_traits
 * \{
 */

// ----------------------------------------------------------------------------
// value_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `value_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::input_iterator.
 */
template <std::input_iterator it_t>
struct value_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_value_t<it_t>;
};

// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `reference` type of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::input_iterator.
 */
template <std::input_iterator it_t>
struct reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_reference_t<it_t>;
};

// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `rvalue_reference` type of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::input_iterator.
 */
template <std::input_iterator it_t>
struct rvalue_reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_rvalue_reference_t<it_t>;
};

// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

// only defined for ranges

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `difference_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::weakly_incrementable.
 */
template <std::weakly_incrementable it_t>
struct difference_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::iter_difference_t<it_t>;
};

// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `size_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type you wish to query; must model std::weakly_incrementable.
 */
template <std::weakly_incrementable it_t>
struct size_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::make_unsigned_t<std::iter_difference_t<it_t>>;
};

// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// iterator_tag
// ----------------------------------------------------------------------------

/*!\brief Exposes the
 * [iterator_category](https://en.cppreference.com/w/cpp/iterator/iterator_tags) from the modelled concept.
 * \implements seqan3::transformation_trait
 * \tparam it_t The type to operate on.
 * \see seqan3::iterator_tag_t
 *
 * \attention
 * If [std::iterator_traits<it_t>::iterator_category](https://en.cppreference.com/w/cpp/iterator/iterator_traits)
 * is defined for a type `it_t`, this transformation trait acts as an alias for it.
 * If it is not defined and no concepts are modelled, iterator_tag<it_t>::type is not defined.
 */
template <typename it_t>
struct iterator_tag
{
SEQAN3_DOXYGEN_ONLY(
    //!\brief The [iterator_category](https://en.cppreference.com/w/cpp/iterator/iterator_tags).
    using type = iterator_category;
)
};

//!\cond
template <typename it_t>
    requires requires { typename std::iterator_traits<it_t>::iterator_category; }
struct iterator_tag<it_t>
{
    using type = typename std::iterator_traits<it_t>::iterator_category;
};

template <std::input_iterator it_t>
    requires (!requires { typename std::iterator_traits<it_t>::iterator_category; })
struct iterator_tag<it_t>
{
    using type = std::input_iterator_tag;
};

template <typename it_t>
    requires (!std::input_iterator<it_t>) && std::output_iterator<it_t, std::iter_value_t<it_t>> &&
             (!requires { typename std::iterator_traits<it_t>::iterator_category; })
struct iterator_tag<it_t>
{
    using type = std::output_iterator_tag;
};

template <std::forward_iterator it_t>
    requires (!requires { typename std::iterator_traits<it_t>::iterator_category; })
struct iterator_tag<it_t>
{
    using type = std::forward_iterator_tag;
};

template <std::bidirectional_iterator it_t>
    requires (!requires { typename std::iterator_traits<it_t>::iterator_category; })
struct iterator_tag<it_t>
{
    using type = std::bidirectional_iterator_tag;
};

template <std::random_access_iterator it_t>
    requires (!requires { typename std::iterator_traits<it_t>::iterator_category; })
struct iterator_tag<it_t>
{
    using type = std::random_access_iterator_tag;
};
//!\endcond

/*!\brief Return the `iterator_category` type of the input type (transformation_trait shortcut).
 * \tparam it_t The type to operate on.
 * \see seqan3::iterator_tag
 */
template <typename it_t>
//!\cond
    requires requires { typename iterator_tag<it_t>::type; }
//!\endcond
using iterator_tag_t = typename iterator_tag<it_t>::type;

//!\}

} // namespace seqan3
