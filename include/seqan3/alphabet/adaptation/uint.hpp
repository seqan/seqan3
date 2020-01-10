// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
#include <seqan3/core/detail/int_types.hpp>

namespace seqan3::detail
{
//!\brief Whether a type is `uint8_t`, `uint16_t` or `uint32_t`.
//!\ingroup adaptation
//!\hideinitializer
template <typename type>
constexpr bool is_uint_adaptation_v = std::same_as<type, uint8_t>  ||
                                      std::same_as<type, uint16_t> ||
                                      std::same_as<type, uint32_t>;
} // namespace seqan3::detail

namespace seqan3::custom
{

/*!\brief Alphabet specific customisations for unsigned integral types.
 * \tparam uint_type Any of `uint8_t`, `uint16_t` and `uint32_t`.
 * \ingroup adaptation
 */
template <typename uint_type>
//!\cond
    requires seqan3::detail::is_uint_adaptation_v<uint_type>
//!\endcond
struct alphabet<uint_type>
{
    //!\brief Return the number of values the uint type can take (e.g. 256 for `uint8_t`).
    static constexpr auto alphabet_size =
        detail::min_viable_uint_t<detail::size_in_values_v<uint_type>>{detail::size_in_values_v<uint_type>};

    /*!\brief Converting uint to char casts to a character type of same size.
     * \param[in] intgr The alphabet letter that you wish to convert to char.
     * \returns The letter's value in the alphabet's rank type (usually `uint`).
     */
    static constexpr auto to_char(uint_type const intgr) noexcept
    {
        if constexpr (std::same_as<uint_type, uint8_t>)
            return static_cast<char>(intgr);
        else if constexpr (std::same_as<uint_type, uint16_t>)
            return static_cast<char16_t>(intgr);
        else
            return static_cast<char32_t>(intgr);
    }

    /*!\brief Converting uint to rank is a no-op (it will just return the value you pass in).
     * \param[in] intgr The alphabet letter that you wish to convert to rank.
     * \returns `intgr`.
     */
    static constexpr uint_type to_rank(uint_type const intgr) noexcept
    {
        return intgr;
    }

    /*!\brief Assign from a character type via implicit or explicit cast.
     * \param[in] chr The `char` value you wish to assign.
     * \param[in,out] intgr The alphabet letter that you wish to assign to.
     * \returns A reference to the alphabet letter you passed in.
     */
    static constexpr uint_type & assign_char_to(decltype(to_char(uint_type{})) const chr, uint_type & intgr) noexcept
    {
        return intgr = chr;
    }

    /*!\brief Assign a rank to to the uint (same as calling `=`).
     * \param[in] intgr2 The `rank` value you wish to assign.
     * \param[in,out] intgr The alphabet letter that you wish to assign to.
     * \returns A reference to the alphabet letter you passed in.
     */
    static constexpr uint_type & assign_rank_to(uint_type const intgr2, uint_type & intgr) noexcept
    {
        return intgr = intgr2;
    }
};

} // namespace seqan3::custom
