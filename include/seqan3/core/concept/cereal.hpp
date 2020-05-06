// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Adaptions of concepts from the Cereal library.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/details/traits.hpp>
#include <cereal/archives/binary.hpp>
#endif

namespace seqan3
{

/*!\interface seqan3::cereal_output_archive <>
 * \brief All output archives of the Cereal library satisfy this.
 * \extends seqan3::cereal_archive
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
SEQAN3_CONCEPT cereal_output_archive = std::is_base_of_v<cereal::detail::OutputArchiveBase, t>;
#else
template <typename t>
SEQAN3_CONCEPT cereal_output_archive = false;
#endif
//!\endcond

/*!\interface seqan3::cereal_input_archive <>
 * \brief All input archives of the Cereal library satisfy this.
 * \extends seqan3::cereal_archive
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
SEQAN3_CONCEPT cereal_input_archive = std::is_base_of_v<cereal::detail::InputArchiveBase, t>;
#else
template <typename t>
SEQAN3_CONCEPT cereal_input_archive = false;
#endif
//!\endcond

/*!\interface seqan3::cereal_archive <>
 * \brief All archives of the Cereal library satisfy this.
 * \ingroup core
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename t>
SEQAN3_CONCEPT cereal_archive = cereal_output_archive<t> || cereal_input_archive<t>;
#else
template <typename t>
SEQAN3_CONCEPT cereal_archive = false;
#endif
//!\endcond

/*!\interface seqan3::cereal_text_archive <>
 * \brief All text archives of the Cereal library satisfy this.
 * \extends seqan3::cereal_archive
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
SEQAN3_CONCEPT cereal_text_archive = std::is_base_of_v<cereal::traits::TextArchive, t>;
#else
template <typename t>
SEQAN3_CONCEPT cereal_text_archive = false;
#endif
//!\endcond

/*!\interface seqan3::cerealisable <>
 * \ingroup core
 * \brief Specifies the requirements for types that are serialisable via Cereal.
 *
 * The `value_t` type satisfy the cerealisable, if `value_t` can be
 * serialised with cereal, i.e. `value_t` has a single serialisation function
 * (`serialize`) or split load/save pair (load and save) either inside or
 * outside of the class.
 *
 * \sa https://uscilab.github.io/cereal/serialization_functions.html
 *
 * ```
 * #include <seqan3/core/concept/cereal.hpp>
 *
 * // fundamental types are serialisable
 * static_assert(seqan3::cerealisable<int>);
 *
 * #include <array>
 * #include <cereal/types/array.hpp> // std::array is now serialisable
 * static_assert(seqan3::cerealisable<std::array<int, 12>>);
 *
 * #include <seqan3/alphabet/nucleotide/dna4.hpp> // dna4 is serialisable
 * static_assert(seqan3::cerealisable<seqan3::dna4>);
 * ```
 *
 * ### Example
 *
 * \include test/snippet/core/cereal_example.cpp
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename value_t,
          typename input_archive_t = cereal::BinaryInputArchive,
          typename output_archive_t = cereal::BinaryOutputArchive>
SEQAN3_CONCEPT cerealisable =
    cereal::traits::is_input_serializable<value_t, input_archive_t>::value &&
    cereal::traits::is_output_serializable<value_t, output_archive_t>::value;
#else
template <typename value_t,
          typename input_archive_t = void,
          typename output_archive_t = void>
SEQAN3_CONCEPT cerealisable = false;
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
