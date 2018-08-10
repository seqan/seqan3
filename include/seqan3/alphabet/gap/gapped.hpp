// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Contains seqan3::gapped.
 */

#pragma once

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/composition/union_composition.hpp>

namespace seqan3
{


/*!\brief Extends a given alphabet with a gap character.
 * \ingroup gap
 * \tparam alphabet_t Type of the letter, e.g. dna4; must satisfy seqan3::alphabet_concept.
 *
 * The gapped alphabet represents the union of a given alphabet and the
 * seqan3::gap alphabet (e.g. the four letter DNA alphabet + a gap character).
 * Note that you cannot assign regular characters, but additional functions for
 * this are available.
 *
 * ```cpp
 * gapped<dna4> gapped_letter{};
 * gapped<dna4> converted_letter{dna4::C};
 * // doesn't work:
 * // gapped<dna4> my_letter{'A'};
 *
 * gapped<dna4>{}.assign_char('C'); // <- this does!
 * gapped<dna4>{}.assign_char('-'); // gap character
 * gapped<dna4>{}.assign_char('K'); // unknown characters map to the default/unknown
 *                                  // character of the given alphabet type (i.e. A of dna4)
 * ```
 *
 * \sa For more details see union_composition, which is the base class and more general than the gapped alphabet.
 */
template <typename alphabet_t>
//!\cond
    requires alphabet_concept<alphabet_t>
//!\endcond
struct gapped : public union_composition<alphabet_t, gap>
{
    using union_composition<alphabet_t, gap>::_value;
    using union_composition<alphabet_t, gap>::value_size;

    using union_composition<alphabet_t, gap>::union_composition;

    using typename union_composition<alphabet_t, gap>::rank_type;
    using typename union_composition<alphabet_t, gap>::char_type;

    //!\copydoc union_composition::assign_char
    constexpr gapped & assign_char(char_type const c) noexcept
    {
        // We can't just use `using union_composition<alphabet_t, gap>::assign_char;` and need to explicitly forward
        // `assign_char`, because otherwise the return type would be `union_composition` and not `gapped`, which is
        // required by the `alphabet_concept`.
        union_composition<alphabet_t, gap>::assign_char(c);
        return *this;
    }

    //!\copydoc union_composition::assign_rank
    constexpr gapped & assign_rank(rank_type const i) /*noexcept*/
    {
        // TODO(marehr): mark function noexcept if assert (within union_composition) is replaced
        // https://github.com/seqan/seqan3/issues/85
        union_composition<alphabet_t, gap>::assign_rank(i);
        return *this;
    }
};

} // namespace seqan3
