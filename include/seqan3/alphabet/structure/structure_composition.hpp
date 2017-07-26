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
 * \ingroup alphabet
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Contains structure alphabet compositions.
 */

#pragma once

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/structure/concept.hpp>

namespace seqan3
{

/*!\brief A seqan3::cartesian_composition that joins a nucleotide alphabet with a structure alphabet.
 * \ingroup alphabet
 * \tparam sequence_alphabet_t Type of the first letter; must satisfy seqan3::nucleotide_concept.
 * \tparam structure_alphabet_t Types of further letters (up to 4); must satisfy seqan3::structure_concept.
 *
 * This composition pairs a nucleotide alphabet with a structure alphabet. The rank values
 * correpsond to numeric values in the size of the composition, while the character values
 * are taken from the sequence alphabet and the structure annotation is taken from the structure
 * alphabet.
 *
 * As with all `seqan3::cartesian_composition` s you may access the individual alphabet letters in
 * regular c++ tuple notation, i.e. `get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * ~~~~~~~~~~~~~~~{.cpp}
 *
 * structure_composition<rna4, dot_bracket3> l{rna4::G, dot_bracket3::PAIR_OPEN};
 * std::cout << int(to_rank(l)) << ' '
 *           << int(to_rank(get<0>(l))) << ' '
 *           << int(to_rank(get<1>(l))) << '\n';
 * // 6 2 1
 *
 * std::cout << to_char(l) << ' '
 *           << to_char(get<0>(l)) << ' '
 *           << to_char(get<1>(l)) << '\n';
 * // G G (
 *
 * // modify via structured bindings and references:
 * auto & [ seq_l, structure_l ] = l;
 * seq_l = rna4::U;
 * std::cout << to_char(l) << '\n';
 * // U
 *
 * ~~~~~~~~~~~~~~~
 *
 * This seqan3::cartesian_composition itself fulfills both seqan3::alphabet_concept and seqan3::structure_concept .
 */

template <typename sequence_alphabet_t, typename structure_alphabet_t>
      requires nucleotide_concept<sequence_alphabet_t> &&
               structure_concept<structure_alphabet_t>
struct structure_composition :
    public cartesian_composition<structure_composition<sequence_alphabet_t, structure_alphabet_t>,
                                 sequence_alphabet_t, structure_alphabet_t>
{
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;
    //!\brief Second template parameter as member type.
    using structure_alphabet_type = structure_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = underlying_char_t<sequence_alphabet_type>;
    //!\brief Equals the char_type of the structure_alphabet_type.
    using str_char_type = underlying_char_t<structure_alphabet_type>;

    /*!\name Write functions
     * \{
     */
    //!\brief Directly assign the sequence character.
    constexpr structure_composition & operator=(sequence_alphabet_type const l) noexcept
    {
        get<0>(*this) = l;
        return *this;
    }

    //!\brief Directly assign the structure character.
    constexpr structure_composition & operator=(structure_alphabet_type const l) noexcept
    {
        get<1>(*this) = l;
        return *this;
    }

    //!\brief Assign from a nucleotide character. This modifies the internal sequence letter.
    constexpr structure_composition & assign_char(char_type const c)
    {
        seqan3::assign_char(get<0>(*this), c);
        return *this;
    }

    //!\brief Assign from a structure character. This modifies the internal structure letter.
    constexpr structure_composition & assign_structure(str_char_type const c)
    {
        seqan3::assign_char(get<1>(*this), c);
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the structure character. This reads the internal structure letter.
    constexpr str_char_type to_structure() const noexcept
    {
        return seqan3::to_char(get<1>(*this));
    }

    //!\brief Return a character. This reads the internal sequence letter.
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }
    //!\}
};

//!\brief Type deduction guide enables usage of structure_composition without specifying template args.
//!\relates structure_composition
template <typename sequence_alphabet_type, typename structure_alphabet_type>
structure_composition(sequence_alphabet_type &&, structure_alphabet_type &&)
    -> structure_composition<std::decay_t<sequence_alphabet_type>, std::decay_t<structure_alphabet_type>>;

} // namespace seqan3

namespace seqan3::detail
{

//!\brief Since seqan3::structure_composition wraps a nucleotide alphabet it is also one.
template <typename sequence_alphabet_type, typename structure_alphabet_type>
struct is_nucleotide<structure_composition<sequence_alphabet_type, structure_alphabet_type>> : public std::true_type
{};

//!\brief Since seqan3::structure_composition wraps a structure alphabet it is also one.
template <typename sequence_alphabet_type, typename structure_alphabet_type>
struct is_structure<structure_composition<sequence_alphabet_type, structure_alphabet_type>> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
static_assert(seqan3::nucleotide_concept<seqan3::structure_composition<seqan3::rna4, seqan3::dot_bracket3>>);
static_assert(seqan3::structure_concept<seqan3::structure_composition<seqan3::rna4, seqan3::dot_bracket3>>);
#endif
