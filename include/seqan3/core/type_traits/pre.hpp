// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various transformation trait base templates and shortcuts.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

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
 * \tparam t The type you wish to query.
 * \see seqan3::value_type_t
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this transformation trait: seqan3::value_type_t
 */
template <typename t>
struct value_type;

/*!\brief Shortcut for seqan3::value_type (transformation_trait shortcut).
 * \tparam t The type you wish to query.
 * \see seqan3::value_type
 */
template <typename t>
using value_type_t = typename value_type<t>::type;

// see specialisation for iterators in core/type_traits/iterator.hpp
// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `reference` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query.
 * \see seqan3::reference_t
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this transformation trait: seqan3::reference_t
 */
template <typename t>
struct reference;

/*!\brief Shortcut for seqan3::reference (transformation_trait shortcut).
 * \tparam t The type you wish to query.
 * \see seqan3::reference
 */
template <typename t>
using reference_t = typename reference<t>::type;

// see specialisation for iterators in core/type_traits/iterator.hpp
// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `rvalue_reference` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query.
 * \see seqan3::rvalue_reference_t
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this transformation trait: seqan3::rvalue_reference_t
 */
template <typename t>
struct rvalue_reference;

/*!\brief Shortcut for seqan3::rvalue_reference (transformation_trait shortcut).
 * \tparam t The type you wish to query.
 * \see seqan3::rvalue_reference
 */
template <typename t>
using rvalue_reference_t = typename rvalue_reference<t>::type;

// see specialisation for iterators in core/type_traits/iterator.hpp
// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `const_reference` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query.
 * \see seqan3::const_reference_t
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this transformation trait: seqan3::const_reference_t
 *
 * \attention This transformation trait is not overloaded for iterators by default, but it is for ranges.
 */
template <typename t>
struct const_reference;

/*!\brief Shortcut for seqan3::const_reference (transformation_trait shortcut).
 * \tparam t The type you wish to query.
 * \see seqan3::const_reference
 */
template <typename t>
using const_reference_t = typename const_reference<t>::type;

// no specialisation for iterators
// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `difference_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query.
 * \see seqan3::difference_type_t
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this transformation trait: seqan3::difference_type_t
 */
template <typename t>
struct difference_type;

/*!\brief Shortcut for seqan3::difference_type (transformation_trait shortcut).
 * \tparam t The type you wish to query.
 * \see seqan3::difference_type
 */
template <typename t>
using difference_type_t = typename difference_type<t>::type;

// see specialisation for iterators in core/type_traits/iterator.hpp
// see specialisation for ranges in core/type_traits/range.hpp

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `size_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query.
 * \see seqan3::size_type_t
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this transformation trait: seqan3::size_type_t
 */
template <typename t>
struct size_type;

/*!\brief Shortcut for seqan3::size_type (transformation_trait shortcut).
 * \tparam t The type you wish to query.
 * \see seqan3::size_type
 */
template <typename t>
using size_type_t = typename size_type<t>::type;

// see specialisation for iterators in core/type_traits/iterator.hpp
// see specialisation for ranges in core/type_traits/range.hpp

//!\}

} // namespace seqan3
