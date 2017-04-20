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

/*!\file alphabet/quality/composition.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains quality alphabet compositions.
 */

#pragma once

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/composition.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace seqan3
{

/*!\brief An alphabet_composition that joins a nucleotide alphabet with a quality alphabet.
 * \ingroup alphabet
 * \tparam sequence_alphabet_t Type of the first letter; must satisfy seqan3::nucleotide_concept.
 * \tparam quality_alphabet_t Types of further letters (up to 4); must satisfy seqan3::quality_concept.
 *
 * This composition pairs a nucleotide alphabet with a quality alphabet. The rank values
 * correpsond to numeric values in the size of the composition, while the character values
 * are taken from the sequence alphabet and the phred values are taken from the quality
 * alphabet.
 *
 * As with all `alphabet_composition` s you may access the individual alphabet letters in
 * regular c++ tuple notation, i.e. `get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * ~~~~~~~~~~~~~~~{.cpp}
 *
 * quality_composition<dna4, illumina18> l{dna4::A, 7};
 * std::cout << int(to_rank(l)) << ' '
 *           << int(to_rank(get<0>(l))) << ' '
 *           << int(to_rank(get<1>(l))) << '\n';
 * // 148 0 7
 *
 * std::cout << to_char(l) << ' '
 *           << to_char(get<0>(l)) << ' '
 *           << to_char(get<1>(l)) << '\n';
 * // A A (
 *
 * std::cout << int(to_phred(l)) << ' '
 * //           << int(to_phred(get<0>(l))) << ' ' // dna4 doesn't have a phred
 *           << int(to_phred(get<1>(l))) << '\n';
 * // 7 7
 *
 * // modify via structured bindings and references:
 * auto & [ seq_l, qual_l ] = l;
 * seq_l = dna4::G;
 * std::cout << to_char(l) << '\n';
 * // G
 *
 * ~~~~~~~~~~~~~~~
 *
 * This alphabet_composition itself fulfills both seqan3::alphabet_concept and seqan3::quality_concept .
 */

template <typename sequence_alphabet_t, typename quality_alphabet_t>
      requires nucleotide_concept<sequence_alphabet_t> &&
               quality_concept<quality_alphabet_t>
struct quality_composition :
    public alphabet_composition<quality_composition<sequence_alphabet_t, quality_alphabet_t>,
                                sequence_alphabet_t, quality_alphabet_t>
{
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;
    //!\brief Second template parameter as member type.
    using quality_alphabet_type = quality_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = underlying_char_t<sequence_alphabet_type>;
    //!\brief Equals the phred_type of the quality_alphabet_type.
    using phred_type = underlying_phred_t<quality_alphabet_type>;

    /*!\name Write functions
     * \{
     */
    //!\brief Directly assign the sequence letter.
    constexpr quality_composition & operator=(sequence_alphabet_type const l) noexcept
    {
        get<0>(*this) = l;
        return *this;
    }

    //!\brief Directly assign the quality letter.
    constexpr quality_composition & operator=(quality_alphabet_type const l) noexcept
    {
        get<1>(*this) = l;
        return *this;
    }

    //!\brief Assign from a character. This modifies the internal sequence letter.
    constexpr quality_composition & assign_char(char_type const c)
    {
        seqan3::assign_char(get<0>(*this), c);
        return *this;
    }

    //!\brief Assign from a phred value. This modifies the internal quality letter.
    constexpr quality_composition & assign_phred(phred_type const c)
    {
        seqan3::assign_phred(get<1>(*this), c);
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the phred value. This reads the internal quality letter.
    constexpr phred_type to_phred() const noexcept
    {
        return seqan3::to_phred(get<1>(*this));
    }

    //!\brief Return a character. This reads the internal sequence letter.
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }
    //!\}
};

//!\brief Type deduction guide enables usage of quality_composition without specifying template args.
//!\relates quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
quality_composition(sequence_alphabet_type &&, quality_alphabet_type &&)
    -> quality_composition<std::decay_t<sequence_alphabet_type>, std::decay_t<quality_alphabet_type>>;

} // namespace seqan3

namespace seqan3::detail
{

//!\brief Since seqan3::quality_composition wraps a nucleotide alphabet it is also one.
template <typename sequence_alphabet_type, typename quality_alphabet_type>
struct is_nucleotide<quality_composition<sequence_alphabet_type, quality_alphabet_type>> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/illumina18.hpp>
static_assert(seqan3::nucleotide_concept<seqan3::quality_composition<seqan3::dna4, seqan3::illumina18>>);
static_assert(seqan3::quality_concept<seqan3::quality_composition<seqan3::dna4, seqan3::illumina18>>);
#endif
