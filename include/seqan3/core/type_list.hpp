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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::type_list and auxiliary metafunctions.
 */

#pragma once

#include <meta/meta.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief Type that contains multiple types, an alias for
 * [meta::list](https://ericniebler.github.io/range-v3/structmeta_1_1list.html).
 * \ingroup core
 */
template <typename ... types>
using type_list = meta::list<types...>;

/*!\brief Type metafunction that extracts the types from a seqan3::type_list and specialises another template with them.
 * \tparam target_type  The type template you wish to specialise.
 * \tparam type_list    The seqan3::type_list with the source types.
 * \ingroup core
 *
 * Enables using the types contained in a seqan3::type_list to specialise another type template. A metafunction
 * shortcut is also defined: seqan3::unpack_type_list_onto_t
 *
 * ### Example
 *
 * ```cpp
 * using tl = type_list<int, char, double>;
 * using t = unpack_type_list_onto_t<std::tuple, tl>;
 * // t is std::tuple<int, char, double>
 * ```
 */

//!\cond
template <template <typename ...> typename target_type, typename type_list_t>
struct unpack_type_list_onto;
//!\endcond

template <template <typename ...> typename target_type, typename ... types>
struct unpack_type_list_onto<target_type, type_list<types...>>
{
    //!\brief The return type: the target type specialised by the unpacked types in the list.
    using type = target_type<types...>;
};

/*!\brief Type metafunction shortcut for seqan3::unpack_type_list_onto.
 * \ingroup core
 * \relates unpack_type_list_onto
 */
template <template <typename ...> typename target_type, typename type_list_t>
using unpack_type_list_onto_t = typename unpack_type_list_onto<target_type, type_list_t>::type;

} // namespace seqan3
