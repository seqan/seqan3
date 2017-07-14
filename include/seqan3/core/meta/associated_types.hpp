// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// ==========================================================================

/*!\file core/meta/associated_types.hpp
 * \brief Meta functions to query type information of associated types.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \ingroup core
 */

#pragma once

#include <meta/meta.hpp>

#include <range/v3/utility/associated_types.hpp>
#include <range/v3/utility/nullptr_v.hpp>

namespace seqan3::detail
{

// default fallback, returns pointer to type
template <typename t, size_t N>
std::remove_cv<t*>
iterator_type_helper(t (*)[N]);

template <typename t>
std::remove_cv<t*>
iterator_type_helper(t **);

template <typename t>
ranges::v3::detail::object_remove_cv<typename t::iterator>
iterator_type_helper(t *);

template <typename t>
ranges::v3::detail::object_remove_cv<typename t::const_iterator>
iterator_type_helper(t const *);

template <typename t>
using iterator_type_ = meta::_t<decltype(detail::iterator_type_helper(ranges::v3::nullptr_v<t>))>;
}  // namespace seqan3::detail

namespace seqan3
{

template <typename t>
struct value_type : ranges::v3::value_type<t>
{};

template <typename t>
using value_type_t = typename value_type<t>::type;

// TODO(rrahn): reference_type

template <typename t>
struct difference_type : meta::_t<ranges::v3::difference_type<t>>
{};

template <typename t>
using difference_type_t = typename difference_type<t>::type;

template <typename t>
struct size_type : meta::_t<ranges::v3::size_type<t>>
{};

template <typename t>
using size_type_t = typename size_type<t>::type;

// TODO(rrahn): iterator_category

template <typename t>
struct iterator_type : meta::defer<detail::iterator_type_, t>
{};

template <typename t>
using iterator_type_t = typename iterator_type<t>::type;

}  // namespace seqan3
