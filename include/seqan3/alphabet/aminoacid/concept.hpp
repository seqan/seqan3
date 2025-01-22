// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Provides seqan3::aminoacid_alphabet.
 */

#pragma once

#include <type_traits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/template_inspection.hpp>

// ============================================================================
// aminoacid_empty_base
// ============================================================================

namespace seqan3
{

/*!\brief This is an empty base class that can be inherited by types that shall model seqan3::aminoacid_alphabet.
 * \ingroup alphabet_aminoacid
 * \see seqan3::enable_aminoacid
 *
 * \stableapi{Since version 3.1.}
 */
struct aminoacid_empty_base
{};

} // namespace seqan3

// ============================================================================
// enable_aminoacid
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void enable_aminoacid(args_t...) = delete;

//!\brief seqan3::detail::customisation_point_object (CPO) definition for seqan3::enable_aminoacid.
//!\ingroup alphabet_aminoacid
template <typename alphabet_t>
struct enable_aminoacid_cpo : public detail::customisation_point_object<enable_aminoacid_cpo<alphabet_t>, 2>
{
    //!\brief CRTP base class seqan3::detail::customisation_point_object.
    using base_t = detail::customisation_point_object<enable_aminoacid_cpo<alphabet_t>, 2>;
    //!\brief Only this class is allowed to import the constructors from #base_t. (CRTP safety idiom)
    using base_t::base_t;

    /*!\brief If `alphabet_type` isn't std::is_nothrow_default_constructible, enable_aminoacid will be called with
     *        std::type_identity instead of a default constructed alphabet.
     */
    template <typename alphabet_type>
    using alphabet_or_type_identity =
        std::conditional_t<std::is_nothrow_default_constructible_v<std::remove_cvref_t<alphabet_type>>
                               && seqan3::is_constexpr_default_constructible_v<std::remove_cvref_t<alphabet_type>>,
                           std::remove_cvref_t<alphabet_type>,
                           std::type_identity<alphabet_type>>;

    /*!\brief CPO overload (check 1 out of 3): explicit customisation via `seqan3::custom::alphabet`
     * \tparam alphabet_type The type of the alphabet. (Needed to defer instantiation for incomplete types.)
     */
    template <typename alphabet_type = alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<2>)(
        /*return*/ std::bool_constant<seqan3::custom::alphabet<alphabet_type>::enable_aminoacid>::value == true /*;*/
    );

    /*!\brief CPO overload (check 2 out of 3): argument dependent lookup (ADL), i.e.
     *        `enable_aminoacid(alphabet_type{})`
     * \tparam alphabet_type The type of the alphabet. (Needed to defer instantiation for incomplete types.)
     *
     * \details
     *
     * If the alphabet_type isn't std::is_nothrow_default_constructible,
     * `enable_aminoacid(std::type_identity<alphabet_type>{})` will be called.
     */
    template <typename alphabet_type = alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<1>)(
        /*return*/ std::bool_constant<enable_aminoacid(alphabet_or_type_identity<alphabet_type>{})>::value == true /*;*/
    );

    /*!\brief CPO overload (check 3 out of 3): `alphabet_type` is derived from seqan3::aminoacid_empty_base
     * \tparam alphabet_type The type of the alphabet. (Needed to defer instantiation for incomplete types.)
     */
    template <typename alphabet_type = alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<0>)(
        /*return*/ std::is_base_of_v<seqan3::aminoacid_empty_base, alphabet_type> == true /*;*/
    );
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\brief A trait that indicates whether a type shall model seqan3::aminoacid_alphabet.
 * \tparam t Type of the argument.
 * \ingroup alphabet_aminoacid
 * \details
 *
 * This is an auxiliary trait that is checked by seqan3::aminoacid_alphabet to verify that a type is an amino acid.
 * This trait should never be read from, instead use seqan3::aminoacid_alphabet.
 * However, user-defined alphabets that want to model seqan3::aminoacid_alphabet need to make sure that it evaluates
 * to `true` for their type.
 *
 * ### Specialisation
 *
 * Do not specialise this trait directly. It acts as a wrapper and looks for three possible implementations
 * (in this order):
 *
 *   1. A `static` member variable `enable_aminoacid` of the class `seqan3::custom::alphabet<t>`.
 *   2. A free function `constexpr bool enable_aminoacid(t) noexcept` in the namespace of your type (or as `friend`).
 *   3. If none of these is found, the default value is defined as:
 *     * `true` if the type inherits from seqan3::aminoacid_empty_base (or seqan3::aminoacid_base),
 *     * `false` otherwise.
 *
 * Implementations of 1. and 2. are required to be marked `constexpr` and the value / return value must be convertible
 * to `bool`.
 * Implementations of 2. are required to be marked `noexcept`. The value passed to functions implementing 2.
 * shall be ignored, it is only used for finding the function via argument-dependent lookup.
 * In case that your type is not seqan3::is_constexpr_default_constructible_v and you wish to provide an implementation
 * for 2., instead overload for `std::type_identity<t>`.
 *
 * To make a type model seqan3::aminoacid_alphabet, it is recommended that you derive from seqan3::aminoacid_base.
 * If that is not possible, choose option 2., and only implement option 1. as a last resort.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/aminoacid/enable_aminoacid.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To change the default behaviour for your own alphabet,
 * follow the above instructions.
 *
 * \experimentalapi{Implementation 2 (free function) is not stable.}
 *
 * \stableapi{Since version 3.1. The name seqan3::enable_aminoacid, Implementation 1,
 *            and Implementation 3 are stable and will not change.}
 */
template <typename t>
inline constexpr bool enable_aminoacid =
    detail::adl_only::enable_aminoacid_cpo<std::remove_cvref_t<t>>::cpo_overload(detail::priority_tag<2>{});

// ============================================================================
// concept
// ============================================================================

/*!\interface seqan3::aminoacid_alphabet <>
 * \extends seqan3::alphabet
 * \brief A concept that indicates whether an alphabet represents amino acids.
 * \ingroup alphabet_aminoacid
 *
 * Since an amino acid alphabet has no specific characteristics (like the complement
 * function for nucleotide alphabets), we distinguish an amino acid alphabet by
 * the seqan3::is_aminoacid type trait.
 *
 * ###Concepts and doxygen
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 *
 * \stableapi{Since version 3.1.}
 */
//!\cond
template <typename type>
concept aminoacid_alphabet = alphabet<type> && enable_aminoacid<type>;
//!\endcond

} // namespace seqan3
