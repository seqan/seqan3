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
 * \brief Adaptions of concepts from the Cereal library.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/details/traits.hpp>
#endif

namespace seqan3
{

/*!\interface seqan3::cereal_output_archive_concept <>
 * \brief All output archives of the Cereal library satisfy this.
 * \extends seqan3::cereal_archive_concept
 * \ingroup core
 *
 * This includes cereal::BinaryOutputArchive, cereal::PortableBinaryOutputArchive, cereal::JSONOutputArchive,
 * and cereal::XMLOutputArchive.
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename t>
concept bool cereal_output_archive_concept = std::is_base_of_v<cereal::detail::OutputArchiveBase, t>;
#else
template <typename t>
concept bool cereal_output_archive_concept = false;
#endif
//!\endcond

/*!\interface seqan3::cereal_input_archive_concept <>
 * \brief All input archives of the Cereal library satisfy this.
 * \extends seqan3::cereal_archive_concept
 * \ingroup core
 *
 * This includes cereal::BinaryInputArchive, cereal::PortableBinaryInputArchive, cereal::JSONInputArchive,
 * and cereal::XMLInputArchive.
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename t>
concept bool cereal_input_archive_concept = std::is_base_of_v<cereal::detail::InputArchiveBase, t>;
#else
template <typename t>
concept bool cereal_input_archive_concept = false;
#endif
//!\endcond

/*!\interface seqan3::cereal_archive_concept <>
 * \brief All archives of the Cereal library satisfy this.
 * \ingroup core
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename t>
concept bool cereal_archive_concept = cereal_output_archive_concept<t> || cereal_input_archive_concept<t>;
#else
template <typename t>
concept bool cereal_archive_concept = false;
#endif
//!\endcond

/*!\interface seqan3::cereal_text_archive_concept <>
 * \brief All text archives of the Cereal library satisfy this.
 * \extends seqan3::cereal_archive_concept
 * \ingroup core
 *
 * This includes cereal::JSONOutputArchive, cereal::XMLOutputArchive, cereal::JSONInputArchive,
 * and cereal::XMLInputArchive.
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename t>
concept bool cereal_text_archive_concept = std::is_base_of_v<cereal::traits::TextArchive, t>;
#else
template <typename t>
concept bool cereal_text_archive_concept = false;
#endif
//!\endcond

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Removes type-mangling that Cereal does with certain types on loading.
 * \details Helpful when defining templatised save/load/serialize functions.
 * \ingroup core
 */
#if SEQAN3_WITH_CEREAL
template <typename type>
using strip_cereal_wrapper_t = typename cereal::traits::strip_minimal<std::decay_t<type>>::type;
#else
template <typename type>
using strip_cereal_wrapper_t = type;
#endif

} // namespace seqan3::detail
