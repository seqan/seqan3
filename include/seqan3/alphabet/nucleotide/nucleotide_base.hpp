// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (chr) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (chr) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::nucleotide_base.
 */

#pragma once

#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that refines seqan3::alphabet_base and is used by the nucleotides.
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 * \tparam char_t       The character type of the alphabet (set this to `void` when defining just a
 *                      seqan3::semi_alphabet_concept).
 */
template <typename derived_type, auto size, typename char_t = char>
class nucleotide_base : public alphabet_base<derived_type, size, char_t>
{
private:
    //!\brief Type of the base class.
    using base_t = alphabet_base<derived_type, size, char_t>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr nucleotide_base() : base_t{} {}
    constexpr nucleotide_base(nucleotide_base const &) = default;
    constexpr nucleotide_base(nucleotide_base &&) = default;
    constexpr nucleotide_base & operator=(nucleotide_base const &) = default;
    constexpr nucleotide_base & operator=(nucleotide_base &&) = default;
    //!\}

    //! Befriend the derived_type so it can instantiate.
    friend derived_type;

public:
    // Import from base:
    using typename base_t::char_type;
    using typename base_t::rank_type;
    using base_t::value_size;
    using base_t::to_rank;

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface.
     */
    //!\{
    static derived_type constexpr A       = assign_char(derived_type{}, 'A');
    static derived_type constexpr B       = assign_char(derived_type{}, 'B');
    static derived_type constexpr C       = assign_char(derived_type{}, 'C');
    static derived_type constexpr D       = assign_char(derived_type{}, 'D');
    static derived_type constexpr G       = assign_char(derived_type{}, 'G');
    static derived_type constexpr H       = assign_char(derived_type{}, 'H');
    static derived_type constexpr K       = assign_char(derived_type{}, 'K');
    static derived_type constexpr M       = assign_char(derived_type{}, 'M');
    static derived_type constexpr N       = assign_char(derived_type{}, 'N');
    static derived_type constexpr R       = assign_char(derived_type{}, 'R');
    static derived_type constexpr S       = assign_char(derived_type{}, 'S');
    static derived_type constexpr T       = assign_char(derived_type{}, 'T');
    static derived_type constexpr U       = assign_char(derived_type{}, 'U');
    static derived_type constexpr V       = assign_char(derived_type{}, 'V');
    static derived_type constexpr W       = assign_char(derived_type{}, 'W');
    static derived_type constexpr Y       = assign_char(derived_type{}, 'Y');
    static derived_type constexpr UNKNOWN = assign_char(derived_type{}, 'N');
    //!\}

    /*!\name Conversion operators
     * \{
     */
//     //!\brief Implicit conversion between dna* and rna* of the same size.
//     //!\tparam other_nucl_type The type to convert to; must satisfy seqan3::nucleotide_concept and have the same \link value_size \endlink.
//     template <typename other_nucl_type>
//     //!\cond
//         requires !std::Same<derived_type, other_nucl_type> && nucleotide_concept<other_nucl_type> &&
//                  value_size == alphabet_size_v<other_nucl_type>
//     //!\endcond
//     constexpr operator other_nucl_type() const noexcept
//     {
//         return other_nucl_type{}.assign_rank(to_rank());
//     }

    //!\brief Explicit conversion to any other nucleotide alphabet (via char representation).
    //!\tparam other_nucl_type The type to convert to; must satisfy seqan3::nucleotide_concept.
    template <typename other_nucl_type>
    //!\cond
        requires !std::Same<derived_type, other_nucl_type> && nucleotide_concept<other_nucl_type>
    //!\endcond
    explicit constexpr operator other_nucl_type() const noexcept
    {
        return detail::convert_through_char_representation<other_nucl_type, derived_type>[to_rank()];
    }
    //!\}

    /*!\name Read functions
     * \{
     */

    /*!\brief Return the complement of the letter.
     *
     * \details
     *
     * See \ref nucleotide for the actual values.
     *
     * Satisfies the seqan3::nucleotide_concept::complement() requirement via the seqan3::complement() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr derived_type complement() const noexcept
    {
        return derived_type::complement_table[to_rank()];
    }
    //!\}
};

} // namespace seqan3
