// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various metafunctions for use on iterators.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <type_traits>

#include <range/v3/utility/iterator_traits.hpp>

#include <seqan3/core/platform.hpp>
#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\addtogroup metafunction
 * \{
 */

// ----------------------------------------------------------------------------
// value_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `value_type` of another type [specialisation for input iterators].
 * \tparam it_t The type you wish to query; must satisfy std::input_Iterator.
 */
template <std::InputIterator it_t>
struct value_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = typename std::iterator_traits<std::remove_reference_t<it_t>>::value_type;
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `reference` of another type [specialisation for input iterators].
 * \tparam it_t The type you wish to query; must satisfy std::input_Iterator.
 */
template <std::InputIterator it_t>
struct reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = typename std::iterator_traits<std::remove_reference_t<it_t>>::reference;
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `rvalue_reference` of another type [specialisation for input iterators].
 * \tparam it_t The type you wish to query; must satisfy std::input_Iterator.
 */
template <std::InputIterator it_t>
struct rvalue_reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = decltype(std::ranges::iter_move(std::declval<it_t &>()));
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

// only defined for ranges

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `difference_type` of another type [specialisation for iterators].
 * \tparam it_t The type you wish to query; must satisfy std::WeaklyIncrementable.
 */
template <std::WeaklyIncrementable it_t>
struct difference_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = typename std::iterator_traits<std::remove_reference_t<it_t>>::difference_type;
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `size_type` of another type [specialisation for iterators].
 * \tparam it_t The type you wish to query; must satisfy std::WeaklyIncrementable.
 */
template <std::WeaklyIncrementable it_t>
struct size_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::make_unsigned_t<difference_type_t<it_t>>;
};

// see specialisation for ranges in core/metafunction/range.hpp

//!\}

} // namespace seqan3
