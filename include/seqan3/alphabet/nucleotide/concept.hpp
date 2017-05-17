// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \ingroup nucleotide
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::nucleotide_concept.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// auxiliary metafunction
// ============================================================================

namespace seqan3::detail
{

//!\brief Metafunction that indicates whether an alphabet is a nucleotide alphabet.
//!\ingroup nucleotide
template <typename type>
struct is_nucleotide : public std::false_type
{};

//!\brief Shortcut for seqan3::detail::is_nucleotide.
//!\ingroup nucleotide
template <typename type>
constexpr bool is_nucleotide_v = is_nucleotide<type>::value;

} // namespace seqan3::detail

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{

//!\brief A concept that indicates whether an alphabet represents nucleotides.
//!\ingroup nucleotide
//!\details Refines the seqan3::alphabet_concept.
template <typename type>
concept bool nucleotide_concept = alphabet_concept<type> && detail::is_nucleotide_v<type>;

} // namespace seqan3

// ============================================================================
// conversion specialisations
// ============================================================================

namespace seqan3::detail
{

/*!\brief Specialisation of seqan3::detail::convert for converting between nucleotide types.
 * \ingroup nucleotide
 * \tparam out_t The type of the output, must satisfy seqan3::nucleotide_concept.
 * \tparam in_t The type of the input, must satisfy seqan3::nucleotide_concept.
 */
template <typename out_t, typename in_t>
    requires nucleotide_concept<out_t> && nucleotide_concept<in_t>
struct convert<out_t, in_t>
{
    //!\brief Implementation of seqan3::convert(); casts the input if possible, otherwise looks in conversion table.
    static constexpr out_t impl(in_t const & in) noexcept
    {
        if constexpr (implicitly_convertible_to_concept<in_t, out_t>)
            return in;
        else if constexpr (explicitly_convertible_to_concept<in_t, out_t>)
            return static_cast<out_t>(in);
        else
            return conversion_table[to_rank(in)];
    }
private:
    //!\brief Whether the type combination needs a conversion table (cannot be cast).
    static constexpr bool can_be_cast = implicitly_convertible_to_concept<in_t, out_t> ||
                                        explicitly_convertible_to_concept<in_t, out_t>;

    //!\brief Construct a conversion table at compile time (converts through character representation).
    static constexpr std::array<out_t, can_be_cast ? 1 : alphabet_size_v<in_t>> conversion_table
    {
        [] () constexpr
        {
            // no table needed if we can cast
            if constexpr (can_be_cast)
            {
                return std::array<out_t, 1>{};
            } else
            {
                std::array<out_t, alphabet_size_v<in_t>> ret{};
                for (std::size_t i = 0; i < alphabet_size_v<in_t>; ++i)
                    assign_char(ret[i], to_char(assign_rank(in_t{}, i)));
                return ret;
            }
        }()
    };
};

} // namespace seqan3::detail
