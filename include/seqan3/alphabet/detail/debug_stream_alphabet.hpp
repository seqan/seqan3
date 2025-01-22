// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
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

/*!\brief The printer used for formatted output of seqan3::alphabet types.
 *
 * Prints the char representation of the given alphabet letter.
 *
 * \tparam alphabet_t The type of the alphabet to be printed.
 * \ingroup alphabet
 */
template <alphabet alphabet_t>
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

/*!\brief The printer used for formatted output of seqan3::mask alphabet.
 *
 * Prints "MASKED" if the letter is masked and "UNMASKED" otherwise.
 *
 * \tparam mask_t The type of the alphabet to be printed. Must be seqan3::mask.
 * \ingroup alphabet_mask
 */
template <std::same_as<mask> mask_t>
struct mask_printer<mask_t>
{
    /*!\brief Print the mask alphabet to the stream
     * \tparam stream_t The type of the stream.
     * \param[in,out] stream The stream to print to.
     * \param[in] arg The mask alphabet letter.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, mask_t const arg) const noexcept
    {
        // seqan3::mask is incomplete at this point, so we cannot use `arg == mask{}`
        stream << (arg == mask_t{} ? "UNMASKED" : "MASKED");
    }
};

} // namespace seqan3
