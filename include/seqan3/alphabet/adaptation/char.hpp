// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides alphabet adaptations for standard char types.
 * \details
 * This file provides function and type trait overloads so that the following types
 * fulfil the seqan3::alphabet:
 *   * `char`
 *   * `char16_t`
 *   * `char32_t`
 *
 * You will likely not use these interfaces directly, they are, however, very helpful
 * for conversions between other alphabets and between other alphabets and characters.
 *
 * \attention
 *   * Note that `signed char` and `unsigned char` are absent from the list, because of
 * their type ambiguity with `int8_t` and `uint8_t`.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/int_types.hpp>

namespace seqan3::detail
{
//!\brief Whether a type is `char`, `char16_t`, `char32_t` or `wchar_t` (type trait).
//!\ingroup adaptation
template <typename type>
constexpr bool is_char_adaptation_v = std::same_as<type, char>     ||
                                      std::same_as<type, char16_t> ||
                                      std::same_as<type, char32_t> ||
                                      std::same_as<type, wchar_t>;
} // namespace seqan3::detail

namespace seqan3::custom
{

/*!\brief Alphabet specific customisations for builtin char types.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \ingroup adaptation
 */
template <typename char_type>
//!\cond
    requires detail::is_char_adaptation_v<char_type>
//!\endcond
struct alphabet<char_type>
{
    //!\brief The number of values the char type can take (e.g. 256 for `char`).
    static constexpr auto alphabet_size =
        detail::min_viable_uint_t<detail::size_in_values_v<char_type>>{detail::size_in_values_v<char_type>};

    /*!\brief Converting char to char is no-op (it will just return the value you pass in).
     * \param[in] chr The alphabet letter that you wish to convert to char (no-op).
     * \returns `chr`.
     */
    static constexpr char_type to_char(char_type const chr) noexcept
    {
        return chr;
    }

    /*!\brief Convert char to rank by casting to an unsigned integral type of same size.
     * \param[in] chr The alphabet letter that you wish to convert to rank.
     * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
     */
    static constexpr auto to_rank(char_type const chr) noexcept
    {
        return static_cast<detail::min_viable_uint_t<alphabet_size - 1>>(chr);
    }

    /*!\brief Assign a char to the char type (same as calling `=`).
     * \param[in] chr2 The `char` value you wish to assign.
     * \param[in,out] chr The alphabet letter that you wish to assign to.
     * \returns A reference to the alphabet letter you passed in.
     */
    static constexpr char_type & assign_char_to(char_type const chr2, char_type & chr) noexcept
    {
        return chr = chr2;
    }

    /*!\brief Assigning a rank to a char is the same as assigning it a numeric value.
     * \param[in] rank The `rank` value you wish to assign.
     * \param[in,out] chr The alphabet letter that you wish to assign to.
     * \returns A reference to the alphabet letter you passed in.
     */
    static constexpr char_type & assign_rank_to(decltype(alphabet::to_rank(char_type{})) const rank,
                                                char_type & chr) noexcept
    {
        return chr = rank;
    }
};

} // namespace seqan3::custom
