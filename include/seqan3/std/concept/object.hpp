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
 * \brief Adaptions of core language concepts.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concept/comparison.hpp>
#include <seqan3/std/concept/core_language.hpp>

namespace seqan3
{

/*!\addtogroup concept
 * \{
 */

/*!\interface   seqan3::object_concept <>
 * \brief       Resolves to `std::is_object_v<t>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_object
 */
//!\cond
template <typename t>
concept bool object_concept =                       std::is_object_v<t>;
//!\endcond

/*!\interface   seqan3::destructible_concept <>
 * \brief       Resolves to `std::is_destructible_v<t>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_destructible
 */
/*!\interface   seqan3::trivially_destructible_concept <>
 * \extends     seqan3::destructible_concept
 * \brief       Resolves to `std::is_trivially_destructible_v<t>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_destructible
 */
/*!\interface   seqan3::nothrow_destructible_concept <>
 * \extends     seqan3::destructible_concept
 * \brief       Resolves to `std::is_nothrow_destructible_v<t>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_destructible
 */
//!\cond
template <typename t>
concept bool destructible_concept =                 std::is_destructible_v<t>;
template <typename t>
concept bool trivially_destructible_concept =       destructible_concept<t> &&
                                                    std::is_trivially_destructible_v<t>;
template <typename t>
concept bool nothrow_destructible_concept =         destructible_concept<t> &&
                                                    std::is_nothrow_destructible_v<t>;
//!\endcond

/*!\interface   seqan3::constructible_concept <>
 * \brief       Resolves to `std::is_constructible_v<t, arg_ts...>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_constructible
 */
/*!\interface   seqan3::trivially_constructible_concept <>
 * \extends     seqan3::constructible_concept
 * \brief       Resolves to `std::is_trivially_constructible_v<t, arg_ts...>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_constructible
 */
/*!\interface   seqan3::nothrow_constructible_concept <>
 * \extends     seqan3::constructible_concept
 * \brief       Resolves to `std::is_nothrow_constructible_v<t, arg_ts...>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_constructible
 */
//!\cond
template <typename t, typename ... arg_ts>
concept bool constructible_concept =                 std::is_constructible_v<t, arg_ts...>;
template <typename t, typename ... arg_ts>
concept bool trivially_constructible_concept =       constructible_concept<t, arg_ts...> &&
                                                     std::is_trivially_constructible_v<t, arg_ts...>;
template <typename t, typename ... arg_ts>
concept bool nothrow_constructible_concept =         constructible_concept<t, arg_ts...> &&
                                                     std::is_nothrow_constructible_v<t, arg_ts...>;
//!\endcond

/*!\interface   seqan3::default_constructible_concept <>
 * \extends     seqan3::constructible_concept
 * \brief       Resolves to `std::is_default_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_default_constructible
 */
/*!\interface   seqan3::trivially_default_constructible_concept <>
 * \extends     seqan3::default_constructible_concept
 * \brief       Resolves to `std::is_trivially_default_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_default_constructible
 */
/*!\interface   seqan3::nothrow_default_constructible_concept <>
 * \extends     seqan3::default_constructible_concept
 * \brief       Resolves to `std::is_nothrow_default_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_default_constructible
 */
//!\cond
template <typename t>
concept bool default_constructible_concept =            constructible_concept<t> &&
                                                        std::is_default_constructible_v<t>;
template <typename t>
concept bool trivially_default_constructible_concept =  default_constructible_concept<t> &&
                                                        std::is_trivially_default_constructible_v<t>;
template <typename t>
concept bool nothrow_default_constructible_concept =    default_constructible_concept<t> &&
                                                        std::is_nothrow_default_constructible_v<t>;
//!\endcond

/*!\interface   seqan3::move_constructible_concept <>
 * \extends     seqan3::constructible_concept
 * \brief       Resolves to `std::is_move_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_move_constructible
 */
/*!\interface   seqan3::trivially_move_constructible_concept <>
 * \extends     seqan3::move_constructible_concept
 * \brief       Resolves to `std::is_trivially_move_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_move_constructible
 */
/*!\interface   seqan3::nothrow_move_constructible_concept <>
 * \extends     seqan3::move_constructible_concept
 * \brief       Resolves to `std::is_nothrow_move_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_move_constructible
 */
//!\cond
template <typename t>
concept bool move_constructible_concept =               constructible_concept<t, t &&> &&
                                                        std::is_move_constructible_v<t>;
template <typename t>
concept bool trivially_move_constructible_concept =     move_constructible_concept<t> &&
                                                        std::is_trivially_move_constructible_v<t>;
template <typename t>
concept bool nothrow_move_constructible_concept =       move_constructible_concept<t> &&
                                                        std::is_nothrow_move_constructible_v<t>;
//!\endcond

/*!\interface   seqan3::copy_constructible_concept <>
 * \extends     seqan3::move_constructible_concept
 * \brief       Resolves to `std::is_copy_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_copy_constructible
 */
/*!\interface   seqan3::trivially_copy_constructible_concept <>
 * \extends     seqan3::copy_constructible_concept
 * \brief       Resolves to `std::is_trivially_copy_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_copy_constructible
 */
/*!\interface   seqan3::nothrow_copy_constructible_concept <>
 * \extends     seqan3::copy_constructible_concept
 * \brief       Resolves to `std::is_nothrow_copy_constructible_v<t>`
 * \sa          http://en.cppreference.com/w/cpp/types/is_copy_constructible
 */
//!\cond
template <typename t>
concept bool copy_constructible_concept =               move_constructible_concept<t> &&
                                                        std::is_copy_constructible_v<t>;
template <typename t>
concept bool trivially_copy_constructible_concept =     copy_constructible_concept<t> &&
                                                        std::is_trivially_copy_constructible_v<t>;
template <typename t>
concept bool nothrow_copy_constructible_concept =       copy_constructible_concept<t> &&
                                                        std::is_nothrow_copy_constructible_v<t>;
//!\endcond

/*!\interface   seqan3::movable_concept
 * \brief       Subsumes seqan3::object_concept, seqan3::move_constructible_concept, seqan3::swappable_concept and
 *              requires that the type be seqan3::assignable_concept from a value of itself.
 * \extends     seqan3::object_concept
 * \extends     seqan3::move_constructible_concept
 * \extends     seqan3::assignable_concept
 * \extends     seqan3::swappable_concept
 *
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Movable
 */
//!\cond
template <typename t>
concept bool movable_concept =                      object_concept<t> &&
                                                    move_constructible_concept<t> &&
                                                    assignable_concept<t &, t> &&
                                                    swappable_concept<t>;
//!\endcond

/*!\interface   seqan3::copyable_concept
 * \brief       Subsumes seqan3::movable_concept, seqan3::copy_constructible_concept, and
 *              requires that the type be seqan3::assignable_concept from a `const &` of itself.
 * \extends     seqan3::movable_concept
 * \extends     seqan3::copy_constructible_concept
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Copyable
 */
/*!\interface   seqan3::trivially_copyable_concept
 * \brief       Resolves to std::is_trivially_copyable_v<t>.
 * \extends     seqan3::copyable_concept
 * \sa          http://en.cppreference.com/w/cpp/types/is_trivially_copyable
 */
//!\cond
template <typename t>
concept bool copyable_concept =                     movable_concept<t> &&
                                                    copy_constructible_concept<t> &&
                                                    assignable_concept<t &, t const &>;
template <typename t>
concept bool trivially_copyable_concept =           copyable_concept<t> &&
                                                    std::is_trivially_copyable_v<t>;
//!\endcond

/*!\interface   seqan3::trivial_concept
 * \brief       A type that satisfies seqan3::trivially_copyable_concept and seqan3::trivially_destructible_concept.
 * \extends     seqan3::trivially_copyable_concept
 * \extends     seqan3::trivially_destructible_concept
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Copyable
 */
//!\cond
template <typename t>
concept bool trivial_concept =                      trivially_copyable_concept<t> &&
                                                    trivially_destructible_concept<t>;
//!\endcond

/*!\interface   seqan3::standard_layout_concept
 * \brief       Resolves to std::is_standard_layout_v<t>.
 * \sa          http://en.cppreference.com/w/cpp/concept/StandardLayoutType
 */
//!\cond
template <typename t>
concept bool standard_layout_concept =              std::is_standard_layout_v<t>;
//!\endcond

/*!\interface   seqan3::semi_regular_concept
 * \brief       Subsumes seqan3::copyable_concept and seqan3::default_constructible_concept.
 * \extends     seqan3::copyable_concept
 * \extends     seqan3::default_constructible_concept
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/SemiRegular
 */
//!\cond
template <typename t>
concept bool semi_regular_concept =                 copyable_concept<t> &&
                                                    default_constructible_concept<t>;
//!\endcond

/*!\interface   seqan3::regular_concept
 * \brief       Subsumes seqan3::semi_regular_concept and seqan3::equality_comparable_concept.
 * \extends     seqan3::semi_regular_concept
 * \extends     seqan3::equality_comparable_concept
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Regular
 */
//!\cond
template <typename t>
concept bool regular_concept =                      semi_regular_concept<t> &&
                                                    equality_comparable_concept<t>;
//!\endcond
//!\}
}  // namespace seqan3
