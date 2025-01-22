// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides alphabet adaptations for standard uint types.
 * \details
 * This file provides function and type trait overloads so that the following types
 * fulfil the seqan3::alphabet:
 *   * `uint8_t`
 *   * `uint16_t`
 *   * `uint32_t`
 *
 * You will likely not use these interfaces directly, they are, however, very helpful
 * for conversions between other alphabets and between other alphabets and characters.
 *
 * \attention
 * Note that `uint64_t` is absent from the list, because there is no corresponding
 * character type.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/detail/integer_traits.hpp>

namespace seqan3::detail
{
//!\brief Whether a type is `uint8_t`, `uint16_t` or `uint32_t`.
//!\ingroup alphabet_adaptation
//!\hideinitializer
template <typename type>
constexpr bool is_uint_adaptation_v =
    std::same_as<type, uint8_t> || std::same_as<type, uint16_t> || std::same_as<type, uint32_t>;
} // namespace seqan3::detail

namespace seqan3::custom
{

/*!\brief Alphabet specific customisations for unsigned integral types.
 * \tparam uint_type Any of `uint8_t`, `uint16_t` and `uint32_t`.
 * \ingroup alphabet_adaptation
 *
 * \stableapi{Since version 3.1.}
 */
template <typename uint_type>
    requires seqan3::detail::is_uint_adaptation_v<uint_type>
struct alphabet<uint_type>
{
    /*!\brief Return the number of values the uint type can take (e.g. 256 for `uint8_t`).
     *
     * \stableapi{Since version 3.1.}
     */
    static constexpr auto alphabet_size =
        detail::min_viable_uint_t<detail::size_in_values_v<uint_type>>{detail::size_in_values_v<uint_type>};

    /*!\brief Converting uint to char casts to a character type of same size.
     * \param[in] uint_v The alphabet letter that you wish to convert to char.
     * \returns The letter's value in the alphabet's rank type (usually `uint`).
     *
     * \stableapi{Since version 3.1.}
     */
    static constexpr auto to_char(uint_type const uint_v) noexcept
    {
        if constexpr (std::same_as<uint_type, uint8_t>)
            return static_cast<char>(uint_v);
        else if constexpr (std::same_as<uint_type, uint16_t>)
            return static_cast<char16_t>(uint_v);
        else
            return static_cast<char32_t>(uint_v);
    }

    /*!\brief Converting uint to rank is a no-op (it will just return the value you pass in).
     * \param[in] uint_v The alphabet letter that you wish to convert to rank.
     * \returns `uint_v`.
     *
     * \stableapi{Since version 3.1.}
     */
    static constexpr uint_type to_rank(uint_type const uint_v) noexcept
    {
        return uint_v;
    }

    /*!\brief Assign from a character type via implicit or explicit cast.
     * \param[in] chr The `char` value you wish to assign.
     * \param[in,out] uint_v The alphabet letter that you wish to assign to.
     * \returns A reference to the alphabet letter you passed in.
     *
     * \stableapi{Since version 3.1.}
     */
    static constexpr uint_type & assign_char_to(decltype(to_char(uint_type{})) const chr, uint_type & uint_v) noexcept
    {
        return uint_v = chr;
    }

    /*!\brief Assign a rank to the uint (same as calling `=`).
     * \param[in] uint2_v The `rank` value you wish to assign.
     * \param[in,out] uint_v The alphabet letter that you wish to assign to.
     * \returns A reference to the alphabet letter you passed in.
     *
     * \stableapi{Since version 3.1.}
     */
    static constexpr uint_type & assign_rank_to(uint_type const uint2_v, uint_type & uint_v) noexcept
    {
        return uint_v = uint2_v;
    }
};

} // namespace seqan3::custom
