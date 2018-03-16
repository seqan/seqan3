// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
 * \brief Provides various metafunctions base templates and shortcuts.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\addtogroup metafunction
 * \{
 */

// ----------------------------------------------------------------------------
// value_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `value_type` of another type [metafunction declaration].
 * \tparam t The type you wish to query.
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this metafunction: seqan3::value_type_t
 */
template <typename t>
struct value_type;

/*!\brief Type metafunction shortcut for seqan3::value_type.
 * \tparam t The type you wish to query.
 */
template <typename t>
using value_type_t = typename value_type<t>::type;

// see specialisation for iterators in core/metafunction/iterator.hpp
// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `reference` of another type [metafunction declaration].
 * \tparam t The type you wish to query.
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this metafunction: seqan3::reference_t
 */
template <typename t>
struct reference;

/*!\brief Type metafunction shortcut for seqan3::reference.
 * \tparam t The type you wish to query.
 */
template <typename t>
using reference_t = typename reference<t>::type;

// see specialisation for iterators in core/metafunction/iterator.hpp
// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `rvalue_reference` of another type [metafunction declaration].
 * \tparam t The type you wish to query.
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this metafunction: seqan3::rvalue_reference_t
 */
template <typename t>
struct rvalue_reference;

/*!\brief Type metafunction shortcut for seqan3::rvalue_reference.
 * \tparam t The type you wish to query.
 */
template <typename t>
using rvalue_reference_t = typename rvalue_reference<t>::type;

// see specialisation for iterators in core/metafunction/iterator.hpp
// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `const_reference` of another type [metafunction declaration].
 * \tparam t The type you wish to query.
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this metafunction: seqan3::const_reference_t
 *
 * \attention This metafunction is not overloaded for iterators by default, but it is for ranges.
 */
template <typename t>
struct const_reference;

/*!\brief Type metafunction shortcut for seqan3::const_reference.
 * \tparam t The type you wish to query.
 */
template <typename t>
using const_reference_t = typename const_reference<t>::type;

// no specialisation for iterators
// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `difference_type` of another type [metafunction declaration].
 * \tparam t The type you wish to query.
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this metafunction: seqan3::difference_type_t
 */
template <typename t>
struct difference_type;

/*!\brief Type metafunction shortcut for seqan3::difference_type.
 * \tparam t The type you wish to query.
 */
template <typename t>
using difference_type_t = typename difference_type<t>::type;

// see specialisation for iterators in core/metafunction/iterator.hpp
// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `size_type` of another type [metafunction declaration].
 * \tparam t The type you wish to query.
 *
 * \details
 *
 * This is a pure declaration, you need to create a *definition* for concrete types
 * or specialized or constrained templates.
 *
 * There is a shortcut for this metafunction: seqan3::size_type_t
 */
template <typename t>
struct size_type;

/*!\brief Type metafunction shortcut for seqan3::size_type.
 * \tparam t The type you wish to query.
 */
template <typename t>
using size_type_t = typename size_type<t>::type;

// see specialisation for iterators in core/metafunction/iterator.hpp
// see specialisation for ranges in core/metafunction/range.hpp

//!\}

} // namespace seqan3
