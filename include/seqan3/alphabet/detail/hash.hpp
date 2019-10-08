// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides overloads for std::hash.
 */

#pragma once

#include <functional>

#include <seqan3/alphabet/concept.hpp>

namespace std
{
/*!\brief Struct for hashing a character.
 * \ingroup alphabet
 * \tparam alphabet_t The type of character to hash; Must model seqan3::semialphabet.
 */
template <typename alphabet_t>
    //!\cond
    requires seqan3::semialphabet<alphabet_t>
    //!\endcond
struct hash<alphabet_t>
{
    /*!\brief Compute the hash for a character.
     * \param[in] character The character to process. Must model seqan3::semialphabet.
     *
     * \returns size_t.
     * \sa seqan3::to_rank.
     */
    size_t operator()(alphabet_t const character) const noexcept
    {
        return seqan3::to_rank(character);
    }
};

/*!\brief Struct for hashing a range of characters.
 * \ingroup alphabet
 * \tparam urng_t The type of the range; Must model std::ranges::input_range and the reference type of the range of the
                  range must model seqan3::semialphabet.
 */
template <ranges::input_range urng_t>
    //!\cond
    requires seqan3::semialphabet<seqan3::reference_t<urng_t>>
    //!\endcond
struct hash<urng_t>
{
    /*!\brief Compute the hash for a range of characters.
     * \param[in] range The input range to process. Must model std::ranges::input_range and the reference type of the
                        range of the range must model seqan3::semialphabet.
     * \returns size_t.
     */
    size_t operator()(urng_t const & range) const noexcept
    {
        using alphabet_t = seqan3::value_type_t<urng_t>;
        size_t result{0};
        hash<alphabet_t> h{};
        for (auto const character : range)
        {
            result *= seqan3::alphabet_size<alphabet_t>;
            result += h(character);
        }
        return result;
    }
};
} // namespace std
