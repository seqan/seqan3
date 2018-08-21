// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Extends a given alphabet with the mask alphabet.
 */

#pragma once

#include <seqan3/alphabet/mask/all.hpp>
#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <locale>

namespace seqan3
{
/*!\brief Implementation of a masked composition, which extends a given alphabet
 * with a mask.
 * \ingroup mask
 * \implements seqan3::alphabet_concept
 * \tparam sequence_alphabet_t Type of the first letter; must satisfy seqan3::semi_alphabet_concept.
 * \tparam mask_t Types of masked letter; must satisfy seqan3::semi_alphabet_concept, defaults to seqan3::mask.
 *
 * \details
 * The masked composition represents a seqan3::cartesian_composition of any given alphabet with the
 * masked alphabet. It allows one to specify which portions of a sequence should be masked,
 * without losing additional information by replacing the sequence directly.
 *
 * \snippet test/snippet/alphabet/mask/masked.cpp general
 */
 template <typename sequence_alphabet_t, typename mask_t = mask>
//!\cond
    requires alphabet_concept<sequence_alphabet_t>
//!\endcond
class masked : public cartesian_composition<masked<sequence_alphabet_t, mask_t>, sequence_alphabet_t, mask_t>
{
private:
    //!\brief The base type.
    using base_type = cartesian_composition<masked<sequence_alphabet_t, mask_t>, sequence_alphabet_t, mask_t>;

public:
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = underlying_char_t<sequence_alphabet_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    masked() = default;
    constexpr masked(masked const &) = default;
    constexpr masked(masked &&) = default;
    constexpr masked & operator =(masked const &) = default;
    constexpr masked & operator =(masked &&) = default;
    ~masked() = default;

    using base_type::base_type; // Inherit non-default constructors

    using base_type::operator=; // Inherit non-default assignment operators

    //!\copydoc cartesian_composition::cartesian_composition(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked(component_type const alph) {} ))
    //!\copydoc cartesian_composition::cartesian_composition(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked(indirect_component_type const alph) {} ))
    //!\copydoc cartesian_composition::operator=(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked & operator=(component_type const alph) {} ))
    //!\copydoc cartesian_composition::operator=(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked & operator=(indirect_component_type const alph) {} ))
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    constexpr masked & assign_char(char_type const c)
    {
        seqan3::assign_char(get<0>(*this), c);
        // TODO: The check of if c is lowercase can be optimized using a lookup table.
        seqan3::assign_rank(get<1>(*this), (c >= 'a' && c <= 'z'));
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return a character.
    constexpr char_type to_char() const noexcept
    {
        if (seqan3::to_rank(get<1>(*this)))
        {
            return std::tolower(seqan3::to_char(get<0>(*this)), std::locale("C"));
        }
        else
        {
            return seqan3::to_char(get<0>(*this));
        }
    }
};

//!\brief Type deduction guide enables usage of masked without specifying template args.
//!\relates masked
template <typename sequence_alphabet_type>
masked(sequence_alphabet_type &&, mask const &)
    -> masked<std::decay_t<sequence_alphabet_type>>;
} //namespace seqan3
