// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\cond DEV
 * \file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::detail::convert_through_char_representation.
 * \endcond
 */

#pragma once

#include <algorithm>
#include <array>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>

// ============================================================================
// conversion to/from char/rank types
// ============================================================================

namespace seqan3::detail
{

/*!\brief A precomputed conversion table for two alphabets based on their char representations.
 * \ingroup alphabet
 * \tparam out_t The type of the output, must satisfy seqan3::Alphabet.
 * \tparam in_t The type of the input, must satisfy seqan3::Alphabet.
 * \hideinitializer
 */
template <Alphabet out_t, Alphabet in_t>
constexpr std::array<out_t, alphabet_size<in_t>> convert_through_char_representation
{
    [] () constexpr
    {
        std::array<out_t, alphabet_size<in_t>> ret{};
        for (typename in_t::rank_type i = 0; i < alphabet_size<in_t>; ++i)
            assign_char_to(to_char(assign_rank_to(i, in_t{})), ret[i]);
        return ret;
    }()
};

} // namespace seqan3::detail
