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
 * \brief Contains the composition of nucleotide with structure alphabets.
 */

#pragma once

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/structure/rna_structure_concept.hpp>

namespace seqan3
{

/*!\brief A seqan3::cartesian_composition that joins a nucleotide alphabet with an RNA structure alphabet.
 * \ingroup structure
 * \tparam sequence_alphabet_t Type of the first letter; must satisfy seqan3::nucleotide_concept.
 * \tparam structure_alphabet_t Types of further letters; must satisfy seqan3::rna_structure_concept.
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
 * structured_rna<rna4, dot_bracket3> l{rna4::G, dot_bracket3::PAIR_OPEN};
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
 * This seqan3::cartesian_composition itself fulfills both seqan3::nucleotide_concept and seqan3::rna_structure_concept.
 */

template <typename sequence_alphabet_t, typename structure_alphabet_t>
    requires nucleotide_concept<sequence_alphabet_t> && rna_structure_concept<structure_alphabet_t>
struct structured_rna :
    public cartesian_composition<structured_rna<sequence_alphabet_t, structure_alphabet_t>,
                                 sequence_alphabet_t, structure_alphabet_t>
{
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;
    //!\brief Second template parameter as member type.
    using structure_alphabet_type = structure_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = underlying_char_t<sequence_alphabet_type>;

    //!\name Write functions
    //!\{

    //!\brief Directly assign the sequence character.
    constexpr structured_rna & operator=(sequence_alphabet_type const l) noexcept
    {
        get<0>(*this) = l;
        return *this;
    }

    //!\brief Directly assign the structure character.
    constexpr structured_rna & operator=(structure_alphabet_type const l) noexcept
    {
        get<1>(*this) = l;
        return *this;
    }

    //!\brief Assign from a nucleotide character. This modifies the internal sequence letter.
    constexpr structured_rna & assign_char(char_type const c)
    {
        seqan3::assign_char(get<0>(*this), c);
        return *this;
    }
    //!\}

    //!\name Read functions
    //!\{

    //!\brief Return a character. This reads the internal sequence letter.
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }

    /*!\brief Return a structured_rna where the sequence letter is converted to its complement.
     * \details
     * See \ref nucleotide for the actual values.
     * Satisfies the seqan3::nucleotide_concept::complement() requirement via the seqan3::complement() wrapper.
     * The structure letter is not modified.
     * \par Complexity
     * Constant.
     * \par Exceptions
     * Guaranteed not to throw.
     */
    constexpr structured_rna complement() const noexcept
    {
        return structured_rna{get<0>(*this).complement(), get<1>(*this)};
    }
    //!\}

    //!\name RNA structure properties
    //!\{

    /*!\brief Check whether the character represents a rightward interaction in an RNA structure.
     * \returns True if the letter represents a rightward interaction, False otherwise.
     */
    constexpr bool is_pair_open() const noexcept
    {
        return get<1>(*this).is_pair_open();
    };

    /*!\brief Check whether the character represents a leftward interaction in an RNA structure.
     * \returns True if the letter represents a leftward interaction, False otherwise.
     */
    constexpr bool is_pair_close() const noexcept
    {
        return get<1>(*this).is_pair_close();
    };

    /*!\brief Check whether the character represents an unpaired position in an RNA structure.
     * \returns True if the letter represents an unpaired site, False otherwise.
     */
    constexpr bool is_unpaired() const noexcept
    {
        return get<1>(*this).is_unpaired();
    };

    //!\brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions.
    static constexpr bool pseudoknot_support{structure_alphabet_t::pseudoknot_support};
    //!\}
};

//!\brief Type deduction guide enables usage of structured_rna without specifying template args.
//!\relates structured_rna
template <typename sequence_alphabet_type, typename structure_alphabet_type>
structured_rna(sequence_alphabet_type &&, structure_alphabet_type &&)
    -> structured_rna<std::decay_t<sequence_alphabet_type>, std::decay_t<structure_alphabet_type>>;

} // namespace seqan3
