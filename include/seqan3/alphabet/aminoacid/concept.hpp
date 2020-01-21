// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Provides seqan3::aminoacid_alphabet.
 */

#pragma once

#include <type_traits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

// ============================================================================
// aminoacid_empty_base
// ============================================================================

namespace seqan3
{

/*!\brief This is an empty base class that can be inherited by types that shall model seqan3::aminoacid_alphabet.
 * \ingroup aminoacid
 * \see seqan3::enable_aminoacid
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
template <typename ...args_t>
void enable_aminoacid(args_t ...) = delete;

//!\brief Customisation point dispatcher for seqan3::enable_aminoacid.
struct enable_aminoacid_dispatcher
{
public:
    // explicit customisation
    SEQAN3_CPO_IMPL(2, std::bool_constant<seqan3::custom::alphabet<strip_type_identity_t<t>>::enable_aminoacid>::value)
    // ADL
    SEQAN3_CPO_IMPL(1, std::bool_constant<enable_aminoacid(t{})>::value)
    // default: derived from base class or not (valid for any type)
    SEQAN3_CPO_IMPL(0, (std::is_base_of_v<seqan3::aminoacid_empty_base, strip_type_identity_t<t>>))

    //!\brief Main dispatching function.
    template <typename alph_t>
    static constexpr bool dispatch() noexcept
    {
        if constexpr (std::is_nothrow_default_constructible_v<alph_t> &&
                      seqan3::is_constexpr_default_constructible_v<alph_t>)
        {
            return impl(priority_tag<2>{}, alph_t{});
        }
        else
        {
            return impl(priority_tag<2>{}, std::type_identity<alph_t>{});
        }
    }
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\brief A trait that indicates whether a type shall model seqan3::aminoacid_alphabet.
 * \tparam t Type of the argument.
 * \ingroup aminoacid
 * \details
 *
 * This is an auxiliary trait that is checked by seqan3::aminoacid_alphabet to verify that a type is an amino acid.
 * This trait should never be read from, instead use seqan3::aminoacid_alphabet.
 * However, user-defined alphabets that want to model seqan3::aminoacid_alphabet need to make sure that it evaluates
 * to `true` for their type.
 *
 * ### Specialisation
 *
 * Do not specialise this trait directly. It acts as a wrapper and looks for two possible implementations
 * (in this order):
 *
 *   1. A `static` member variable `enable_aminoacid` of the class `seqan3::custom::alphabet<t>`.
 *   2. A free function `constexpr bool enable_aminoacid(t) noexcept` in the namespace of your type (or as `friend`).
 *
 * If none of these is found, the default value is defined as:
 *
 *   * `true` if the type inherits from seqan3::aminoacid_empty_base (or seqan3::aminoacid_base),
 *   * `false` otherwise.
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
 */
template <typename t>
inline constexpr bool enable_aminoacid = detail::adl_only::enable_aminoacid_dispatcher::dispatch<remove_cvref_t<t>>();

// ============================================================================
// concept
// ============================================================================

/*!\interface seqan3::aminoacid_alphabet <>
 * \extends seqan3::alphabet
 * \brief A concept that indicates whether an alphabet represents amino acids.
 * \ingroup aminoacid
 *
 * Since an amino acid alphabet has no specific characteristics (like the complement
 * function for nucleotide alphabets), we distinguish an amino acid alphabet by
 * the seqan3::is_aminoacid type trait.
 *
 * ###Concepts and doxygen
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT aminoacid_alphabet = alphabet<type> && enable_aminoacid<type>;
//!\endcond

} // namespace seqan3
