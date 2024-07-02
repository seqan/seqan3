// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/io/stream/concept.hpp>

namespace seqan3
{

/*!
 * \brief The printer used for formatted output of seqan3::alphabet types.
 *
 * Prints the char representation of the given alphabet letter.
 *
 * \tparam alphabet_t The type of the alphabet to be printed.
 * \ingroup alphabet
 */
template <typename alphabet_t>
    requires alphabet<alphabet_t>
struct alphabet_printer<alphabet_t>
{
    /*!\brief Print the alphabet to the stream
     * \tparam stream_t The type of the stream.
     * \param[in,out] stream The stream to print to.
     * \param[in] letter The alphabet letter.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, alphabet_t const letter) const noexcept
    {
        stream << to_char(letter);
    }
};

// forward declare seqan3::mask
class mask;

/*!
 * \brief The printer used for formatted output of seqan3::mask alphabet.
 *
 * Prints "MASKED" if the letter is masked and "UNMASKED" otherwise.
 *
 * \tparam alphabet_t The type of the alphabet to be printed.
 * \ingroup alphabet_mask
 */
template <typename alphabet_t>
    requires std::same_as<std::remove_cvref_t<alphabet_t>, mask>
struct mask_printer<alphabet_t>
{
    /*!\brief Print the mask alphabet to the stream
     * \tparam stream_t The type of the stream.
     * \param[in,out] stream The stream to print to.
     * \param[in] letter The mask alphabet letter.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, alphabet_t const letter) const noexcept
    {
        stream << (letter == std::remove_cvref_t<alphabet_t>{} ? "UNMASKED" : "MASKED");
    }
};

} // namespace seqan3
