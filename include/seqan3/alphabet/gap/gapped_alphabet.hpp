// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

/*!\file
 * \ingroup alphabet
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Contains seqan3::gapped_alphabet.
 */

#pragma once

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/union_alphabet.hpp>

namespace seqan3
{


/*!\brief A gapped_alphabet that extends a given alphabet with a gap character.
 * \ingroup alphabet
 * \tparam alphabet_t Type of the letter, e.g. dna4; must satisfy seqan3::alphabet_concept.
 *
 * The gapped_alphabet represents the union of a given alphabet and the
 * seqan3::gap alphabet (e.g. the four letter DNA alphabet + a gap character).
 * Note that you cannot assign regular characters, but additional functions for
 * this are available.
 *
 * ```cpp
 * gapped_alphabet<dna4> gapped_letter{};
 * gapped_alphabet<dna4> converted_letter{dna4::C};
 * // doesn't work:
 * // gapped_alphabet<dna4> my_letter{'A'};
 *
 * gapped_alphabet<dna4>{}.assign_char('C'); // <- this does!
 * gapped_alphabet<dna4>{}.assign_char('-'); // gap character
 * gapped_alphabet<dna4>{}.assign_char('K'); // unknown characters map to the default/unknown
 *                                           // character of the given alphabet type (i.e. A of dna4)
 * ```
 *
 * \sa For more details see union_alphabet, which is the base class and more general than the gapped_alphabet.
 */
template <typename alphabet_t>
    requires alphabet_concept<alphabet_t>
struct gapped_alphabet : public union_alphabet<alphabet_t, gap>
{
    using union_alphabet<alphabet_t, gap>::_value;
    using union_alphabet<alphabet_t, gap>::value_size;

    using union_alphabet<alphabet_t, gap>::union_alphabet;

    using typename union_alphabet<alphabet_t, gap>::rank_type;
    using typename union_alphabet<alphabet_t, gap>::char_type;

    /*!\brief Returns true if it is a gap
     * \details
     * ```cpp
     * gapped_alphabet<dna4> letter = dna4::T;
     *
     * if (!letter.is_gap())
     *     std::cout << "T is NOT a gap character";
     *
     * letter.set_gap();
     * if (letter.is_gap())
     *     std::cout << "Now it is a gap character";
     * ```
     */
    constexpr bool is_gap() const
    {
        return _value == value_size - 1;
    }

    /*!\brief Change it into a gap.
     * \details
     * ```cpp
     * gapped_alphabet<dna4> letter;
     * letter.set_gap();
     *
     * // the same as set_gap()
     * letter = gap::GAP;
     * ```
     */
    constexpr gapped_alphabet set_gap()
    {
        _value = value_size - 1;
        return *this;
    }

    //!\copydoc union_alphabet::assign_rank
    constexpr gapped_alphabet & assign_rank(rank_type const i)
    {
        union_alphabet<alphabet_t, gap>::assign_rank(i);
        return *this;
    }

    //!\copydoc union_alphabet::assign_char
    constexpr gapped_alphabet & assign_char(char_type const c)
    {
        union_alphabet<alphabet_t, gap>::assign_char(c);
        return *this;
    }
};

} // namespace seqan3

#ifndef NDEBUG
#include <seqan3/alphabet/nucleotide/dna4.hpp>
static_assert(seqan3::alphabet_concept<seqan3::gapped_alphabet<seqan3::dna4>>);
#endif
