// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::transformation_trait_or.
 */

#pragma once

#include <type_traits>

#include <meta/meta.hpp>

namespace seqan3::detail
{

/*!\brief This gives a fallback type if *type_t::type* is not defined.
 * \ingroup metafunction
 * \tparam type_t    The type to use if *type_t::type* is defined.
 * \tparam default_t The type to use otherwise.
 *
 * \details
 *
 * Gives *type_t* back if *T::type* is a member type, otherwise *struct{using type = default_t}*.
 *
 * \include test/snippet/core/metafunction/transformation_trait_or.cpp
 *
 * \attention This might get removed if one of our used libraries offers the same
 * functionality.
 *
 * \par Helper types
 *   seqan3::detail::transformation_trait_or_t as a shorthand for *seqan3::detail::transformation_trait_or::type*
 */
template <typename type_t, typename default_t>
using transformation_trait_or = std::conditional_t<meta::is_trait<type_t>::value,  // check if type_t::type exists
                                                   type_t, // if yes, return type_t
                                                   // otherwise return struct{using type = default_t};
                                                   std::enable_if<true, default_t>>;

/*!\brief Helper type of seqan3::detail::transformation_trait_or
 * \ingroup metafunction
 * \relates transformation_trait_or
 */
template <typename type_t, typename default_t>
using transformation_trait_or_t = typename transformation_trait_or<type_t, default_t>::type;

} // namespace seqan3::detail
