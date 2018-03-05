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
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Contains the composition of aminoacid with structure alphabets.
 */

#pragma once

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/structure/dssp9.hpp>

namespace seqan3
{

/*!\brief A seqan3::cartesian_composition that joins an aminoacid alphabet with a protein structure alphabet.
 * \ingroup structure
 * \tparam sequence_alphabet_t Type of the first aminoacid letter; must satisfy seqan3::alphabet_concept.
 * \tparam structure_alphabet_t Types of further structure letters; must satisfy seqan3::alphabet_concept.
 *
 * This composition pairs an aminoacid alphabet with a structure alphabet. The rank values
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
 * structured_aa<aa27, dssp9> l{aa27::W, dssp9::B};
 * std::cout << int(to_rank(l)) << ' '
 *           << int(to_rank(get<0>(l))) << ' '
 *           << int(to_rank(get<1>(l))) << '\n';
 * // 49 22 1
 *
 * std::cout << to_char(l) << ' '
 *           << to_char(get<0>(l)) << ' '
 *           << to_char(get<1>(l)) << '\n';
 * // W W B
 *
 * // modify via structured bindings and references:
 * auto & [ seq_l, structure_l ] = l;
 * seq_l = aa27::V;
 * std::cout << to_char(l) << '\n';
 * // V
 *
 * ~~~~~~~~~~~~~~~
 *
 * This seqan3::cartesian_composition itself fulfills seqan3::alphabet_concept.
 */

template <typename sequence_alphabet_t = aa27, typename structure_alphabet_t = dssp9>
    requires alphabet_concept<sequence_alphabet_t> && alphabet_concept<structure_alphabet_t>
struct structured_aa :
    public cartesian_composition<structured_aa<sequence_alphabet_t, structure_alphabet_t>,
                                 sequence_alphabet_t, structure_alphabet_t>
{
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;
    //!\brief Second template parameter as member type.
    using structure_alphabet_type = structure_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = underlying_char_t<sequence_alphabet_type>;

    /*!\name Write functions
     * \{
     */
    //!\brief Directly assign the sequence character.
    constexpr structured_aa & operator=(sequence_alphabet_type const l) noexcept
    {
        get<0>(*this) = l;
        return *this;
    }

    //!\brief Directly assign the structure character.
    constexpr structured_aa & operator=(structure_alphabet_type const l) noexcept
    {
        get<1>(*this) = l;
        return *this;
    }

    //!\brief Assign from a nucleotide character. This modifies the internal sequence letter.
    constexpr structured_aa & assign_char(char_type const c)
    {
        seqan3::assign_char(get<0>(*this), c);
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return a character. This reads the internal sequence letter.
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }
    //!\}
};

//!\brief Type deduction guide enables usage of structured_aa without specifying template args.
//!\relates structured_aa
template <typename sequence_alphabet_type, typename structure_alphabet_type>
structured_aa(sequence_alphabet_type &&, structure_alphabet_type &&)
    -> structured_aa<std::decay_t<sequence_alphabet_type>, std::decay_t<structure_alphabet_type>>;

} // namespace seqan3
