// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\cond DEV
 * \file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::detail::convert_through_char_representation.
 * \endcond
 */

#pragma once

#include <array>

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// conversion to/from char/rank types
// ============================================================================

namespace seqan3::detail
{

// clang-format off
/*!\brief A precomputed conversion table for two alphabets based on their char representations.
 * \ingroup alphabet
 * \tparam in_t The type of the input, must satisfy seqan3::alphabet.
 * \tparam out_t The type of the output, must satisfy seqan3::alphabet.
 * \hideinitializer
 */
template <alphabet in_t, alphabet out_t>
constexpr std::array<out_t, alphabet_size<in_t>> convert_through_char_representation
{
    []() constexpr {
        std::array<out_t, alphabet_size<in_t>> ret{};

        // for (decltype(alphabet_size<in_t>) i = 0; ...) causes indefinite compilation :(
        for (auto i = decltype(alphabet_size<in_t>){0}; i < alphabet_size<in_t>; ++i)
            assign_char_to(to_char(assign_rank_to(i, in_t{})), ret[i]);

        return ret;
    }()
};
// clang-format on

} // namespace seqan3::detail
