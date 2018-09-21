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
// ============================================================================

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Contains seqan3::phred68legacy quality scores.
 */

#pragma once

#include <algorithm>
#include <cassert>

#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/core/platform.hpp>

// ------------------------------------------------------------------
// phred68legacy
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief Quality type for Solexa and deprecated Illumina formats.
 * \implements seqan3::quality_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \ingroup quality
 *
 * \details
 *
 * The phred68legacy quality alphabet represents the -5-based phred score range
 * [-5..62] mapped to the ASCII range [';' .. '~']. It represents the Solexa and
 * the Illumina [1.0;1.8[ standard.
 *
 * \snippet test/snippet/alphabet/quality/phred68legacy.cpp general
 */
struct phred68legacy
{
    /*!\name Member types
    * \{
    */
    //!\brief The 0-based rank representation of a quality score.
    using rank_type = uint8_t;

    //!\brief The -5-based integer representation of a quality score assignable with =operator.
    using phred_type = int8_t;

    //!\brief The ';'-based character representation of a quality score.
    using char_type = char;
    //!\}

    //!\privatesection
    /*!\name Member variables.
    * \{
    */
    //!\brief The internal 0-based rank value.
    rank_type _value;

    //!\brief The phred score range size for the Solexa standard.
    static constexpr rank_type value_size{68};

    //!\brief The projection offset between char and rank score representation.
    static constexpr char_type offset_char{';'};

    //! projection offsets of a phred quality score
    static constexpr phred_type offset_phred{-5};
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the letter as a character of char_type.
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::to_char() requirement via the seqan3::to_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr char_type to_char() const noexcept
    {
        return static_cast<char_type>(_value + offset_char);
    }

    /*!\brief Return the letter's numeric phred code.
     *
     * \details
     *
     * Satisfies the seqan3::detail::quality_concept::to_phred() requirement via the seqan3::to_phred() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr phred_type to_phred() const noexcept
    {
        return _value + offset_phred;
    }

    /*!\brief Return the letter as a character of char_type.
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::to_char() requirement via the seqan3::to_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr rank_type to_rank() const noexcept
    {
        return _value;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
     /*!\brief Assign from a character.
      *
      * \details
      *
      * Satisfies the seqan3::alphabet_concept::assign_char() requirement via the seqan3::assign_char() wrapper.
      *
      * \par Complexity
      *
      * Constant.
      */
    constexpr phred68legacy & assign_char(char_type const c)
    {
        assert(c >= offset_char && c < offset_char + value_size);
        _value = char_to_value[c];
        return *this;
    }

    /*!\brief Assign from the numeric phred value.
     *
     * \details
     *
     * Satisfies the seqan3::quality_concept::assign_phred() requirement via the seqan3::assign_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     */
    constexpr phred68legacy & assign_phred(phred_type const p)
    {
        // p >= 0 is implicitly always true
        assert((p >= offset_phred) && (p < value_size + offset_phred));
        _value = p - offset_phred;
        return *this;
    }

    /*!\brief Assign from a the numeric rank value.
     *
     * \details
     *
     * Satisfies the seqan3::semi_alphabet_concept::assign_rank() requirement via the seqan3::assign_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     */
    constexpr phred68legacy & assign_rank(rank_type const p)
    {
        // p implicitly always true due to unsignedness of p
        assert(p < value_size);
        _value = p;
        return *this;
    }
    //!\}

    /*!\name Conversion operators
     * \{
     */
    //!\brief Explicit conversion to any other nucleotide alphabet (via char representation).
    //!\tparam other_nucl_type The type to convert to; must satisfy seqan3::quality_concept.
    //  phred scores below 0 will be truncated, i.e. mapped to phred64::phred_type = 0.
    template <quality_concept other_qual_type>
    explicit constexpr operator other_qual_type() const //noexcept
    {
        return detail::convert_through_phred_representation<other_qual_type, std::decay_t<decltype(*this)>>[std::max<phred_type>(0, to_phred())];
    }
    //!\}

    /*!\name Comparison operators.
    * \{
    */
    constexpr bool operator==(const phred68legacy & rhs) const
    {
        return this->_value == rhs._value;
    }

    constexpr bool operator!=(const phred68legacy & rhs) const
    {
        return this->_value != rhs._value;
    }

    constexpr bool operator<(const phred68legacy & rhs) const
    {
        return this->_value < rhs._value;
    }

    constexpr bool operator>(const phred68legacy & rhs) const
    {
        return this->_value > rhs._value;
    }

    constexpr bool operator<=(const phred68legacy & rhs) const
    {
        return this->_value <= rhs._value;
    }

    constexpr bool operator>=(const phred68legacy & rhs) const
    {
        return this->_value >= rhs._value;
    }
    //!\}

protected:
    //!\privatesection
    //!\copydoc seqan3::phred42::char_to_value
    static constexpr std::array<char_type, 256> char_to_value
    {
        [] () constexpr
        {
            std::array<char_type, 256> ret{};
            for (char_type c = ';'; c <= '~'; ++c)
                ret[c] = c - offset_char;
            return ret;
        }()
    };
};

} // namespace seqan3
