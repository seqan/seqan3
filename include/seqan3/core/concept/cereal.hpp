// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

/*!\interface seqan3::CerealOutputArchive <>
 * \brief All output archives of the Cereal library satisfy this.
 * \extends seqan3::CerealArchive
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
SEQAN3_CONCEPT CerealOutputArchive = std::is_base_of_v<cereal::detail::OutputArchiveBase, t>;
#else
template <typename t>
SEQAN3_CONCEPT CerealOutputArchive = false;
#endif
//!\endcond

/*!\interface seqan3::CerealInputArchive <>
 * \brief All input archives of the Cereal library satisfy this.
 * \extends seqan3::CerealArchive
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
SEQAN3_CONCEPT CerealInputArchive = std::is_base_of_v<cereal::detail::InputArchiveBase, t>;
#else
template <typename t>
SEQAN3_CONCEPT CerealInputArchive = false;
#endif
//!\endcond

/*!\interface seqan3::CerealArchive <>
 * \brief All archives of the Cereal library satisfy this.
 * \ingroup core
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename t>
SEQAN3_CONCEPT CerealArchive = CerealOutputArchive<t> || CerealInputArchive<t>;
#else
template <typename t>
SEQAN3_CONCEPT CerealArchive = false;
#endif
//!\endcond

/*!\interface seqan3::CerealTextArchive <>
 * \brief All text archives of the Cereal library satisfy this.
 * \extends seqan3::CerealArchive
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
SEQAN3_CONCEPT CerealTextArchive = std::is_base_of_v<cereal::traits::TextArchive, t>;
#else
template <typename t>
SEQAN3_CONCEPT CerealTextArchive = false;
#endif
//!\endcond

/*!\interface seqan3::Cerealisable <>
 * \ingroup core
 * \brief Specifies the requirements for types that are serialisable via Cereal.
 *
 * The `value_t` type satisfy the Cerealisable, if `value_t` can be
 * serialised with cereal, i.e. `value_t` has a single serialisation function
 * (`serialize`) or split load/save pair (load and save) either inside or
 * outside of the class.
 *
 * \sa https://uscilab.github.io/cereal/serialization_functions.html
 *
 * ```
 * #include <seqan3/core/concept/cereal.hpp>
 * using namespace seqan3;
 *
 * // fundamental types are serialisable
 * static_assert(Cerealisable<int>);
 *
 * #include <array>
 * #include <cereal/types/array.hpp> // std::array is now serialisable
 * static_assert(Cerealisable<std::array<int, 12>>);
 *
 * #include <seqan3/alphabet/nucleotide/dna4.hpp> // dna4 is serialisable
 * static_assert(Cerealisable<dna4>);
 * ```
 *
 * \attention
 * The cereal library is an optional dependency of SeqAn, if it is not found **no types** satisfy this concept.
 */
//!\cond
#if SEQAN3_WITH_CEREAL
template <typename value_t,
          typename input_archive_t = cereal::BinaryInputArchive,
          typename output_archive_t = cereal::BinaryOutputArchive>
SEQAN3_CONCEPT Cerealisable =
    cereal::traits::is_input_serializable<value_t, input_archive_t>::value &&
    cereal::traits::is_output_serializable<value_t, output_archive_t>::value;
#else
template <typename value_t,
          typename input_archive_t = void,
          typename output_archive_t = void>
SEQAN3_CONCEPT Cerealisable = false;
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
