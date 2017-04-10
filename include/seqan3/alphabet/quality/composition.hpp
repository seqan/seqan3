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
// Author: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
// ============================================================================

#pragma once

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/alphabet.hpp>
#include <seqan3/alphabet/composition.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/illumina18.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

/*!\file alphabet/quality/composition.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains quality alphabet compositions.
 */

namespace seqan3
{

/*!\brief An alphabet_composition that joins a regular alphabet with a sequence alphabet.
 * \ingroup alphabet
 * \tparam sequence_alphabet_type Type of the first letter, e.g. dna4; must satisfy alphabet_concept.
 * \tparam quality_alphabet_type Types of further letters (up to 4); must satisfy quality_concept.
 *
 * This composition pairs a regular alphabet with a quality alphabet. The integral values
 * correpsond to numeric values in the size of the composition, while the character values
 * are taken from the sequence alphabet and the phred values are taken from the quality
 * alphabet.
 *
 * As with all alphabet_composition|s you ma access the individual alphabet letters in
 * regular c++ tuple notation, i.e. `std::get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * TODO example
 *
 *
 * This alphabet_composition itself fulfills both the alphabet_concept and the quality_concept.
 */

template <typename sequence_alphabet_type, typename quality_alphabet_type>
      requires alphabet_concept<sequence_alphabet_type> &&
               quality_concept<quality_alphabet_type>
struct quality_composition : public alphabet_composition<sequence_alphabet_type, quality_alphabet_type>
{
//     using char_type = underlying_char_t<sequence_alphabet_type>;
//     using phred_type = underlying_phred_t<quality_alphabet_type>;
//
//     //!
//     // TODO why no constexpr?
// //     quality_composition & from_integral(integral_type const i)
// //     {
// // //         alphabet_composition<sequence_alphabet_type, quality_alphabet_type>::from_integral(i);
// //         return *this;
// //     }
//
//     constexpr quality_composition & from_char(char_type const c)
//     {
//         seqan3::from_char(std::get<0>(*this), c);
//         return *this;
//     }
//
//     constexpr quality_composition & from_phred(phred_type const c)
//     {
//         seqan3::from_phred(std::get<1>(*this), c);
//         return *this;
//     }
//
//     constexpr phred_type to_phred() const
//     {
//         return seqan3::to_phred(std::get<1>(*this));
//     }
//
//     constexpr char_type to_char() const
//     {
//         return seqan3::to_char(std::get<0>(*this));
//     }
};


//!\brief The type of value_size and `alphabet_size_v<alphabet_composition<...>>`
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
struct underlying_char<quality_composition<sequence_alphabet_type, quality_alphabet_type>>
{
    using type = underlying_char_t<sequence_alphabet_type>;
};

//!\brief The type of value_size and `alphabet_size_v<alphabet_composition<...>>`
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
struct underlying_phred<quality_composition<sequence_alphabet_type, quality_alphabet_type>>
{
    using type = underlying_phred_t<quality_alphabet_type>;
};

//!\brief TODO
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
constexpr quality_composition<sequence_alphabet_type, quality_alphabet_type> &
from_integral(quality_composition<sequence_alphabet_type, quality_alphabet_type> & c,
              underlying_integral_t<quality_composition<sequence_alphabet_type, quality_alphabet_type>> const i)
{
    from_integral(static_cast<alphabet_composition<sequence_alphabet_type, quality_alphabet_type>>(c), i);
    return c;
}

//!\brief TODO
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
constexpr quality_composition<sequence_alphabet_type, quality_alphabet_type> &
from_char(quality_composition<sequence_alphabet_type, quality_alphabet_type> & c,
          underlying_char_t<quality_composition<sequence_alphabet_type, quality_alphabet_type>> const i)
{
    from_char(std::get<0>(c), i);
    return c;
}

//!\brief TODO
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
constexpr quality_composition<sequence_alphabet_type, quality_alphabet_type> &
from_phred(quality_composition<sequence_alphabet_type, quality_alphabet_type> & c,
           underlying_phred_t<quality_composition<sequence_alphabet_type, quality_alphabet_type>> const i)
{
    from_phred(std::get<1>(c), i);
    return c;
}

//!\brief TODO
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
constexpr underlying_integral_t<quality_composition<sequence_alphabet_type, quality_alphabet_type>>
to_integral(quality_composition<sequence_alphabet_type, quality_alphabet_type> const & c)
{
    to_integral(static_cast<alphabet_composition<sequence_alphabet_type, quality_alphabet_type> const &>(c));
    return c;
}


//!\brief TODO
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
constexpr underlying_char_t<quality_composition<sequence_alphabet_type, quality_alphabet_type>>
to_char(quality_composition<sequence_alphabet_type, quality_alphabet_type> const & c)
{
    to_char(std::get<0>(c));
}

//!\brief TODO
//!\memberof quality_composition
template <typename sequence_alphabet_type, typename quality_alphabet_type>
constexpr underlying_phred_t<quality_composition<sequence_alphabet_type, quality_alphabet_type>>
to_phred(quality_composition<sequence_alphabet_type, quality_alphabet_type> const & c)
{
    to_phred(std::get<1>(c));
}



//! An alphabet that stores a dna4 letter and an illumina18 letter at each position.
using dna4q = quality_composition<dna4, illumina18>;

//! An alphabet that stores a dna5 letter and an illumina18 letter at each position.
using dna5q = quality_composition<dna5, illumina18>;

// using rna4q = quality_composition<rna4, illumina18>;
// using rna5q = quality_composition<rna5, illumina18>;

} // namespace seqan3

#ifndef NDEBUG
// static_assert(seqan3::alphabet_concept<seqan3::dna4q>);
// static_assert(seqan3::detail::internal_alphabet_concept<seqan3::dna4q>);
// static_assert(seqan3::quality_concept<seqan3::dna4q>);
// static_assert(seqan3::detail::internal_quality_concept<seqan3::dna4q>);
#endif


