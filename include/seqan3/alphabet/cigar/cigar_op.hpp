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
 * \brief Introduces the cigar_op alphabet.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <type_traits>

#include <seqan3/alphabet/detail/alphabet_base.hpp>

// ------------------------------------------------------------------
// cigar_op
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The (extended) cigar operation alphabet of M,D,I,H,N,P,S,X,=.
 * \ingroup cigar
 *
 * \details
 *
 * The CIGAR string can be either basic or extended. The only difference in the
 * extended cigar alphabet is that aligned bases are classified as an actual
 * match ('=') or mismatch ('X'). In contrast, the basic cigar alphabet
 * only indicated the aligned status with an 'M', without further
 * information if the bases are actually equal or not.
 *
 * The main purpose of the seqan3::cigar_op alphabet is to be used in the
 * seqan3::cigar composition, where a cigar operation is paired with a count
 * value.
 *
 * Example usage:
 * \snippet test/snippet/alphabet/cigar/cigar_op.cpp general
 *
 * \note Usually you do not want to manipulate cigar elements and vectors on
 *       your own but convert an alignment to a cigar and back. See
 *       seqan3::get_cigar_vector for how to convert two aligned sequences into
 *       an cigar_vector.
 */
class cigar_op : public alphabet_base<cigar_op, 9, char>
{
private:
	//!\brief The base class.
	using base_t = alphabet_base<cigar_op, 9, char>;

	//!\cond \brief Befriend seqan3::alphabet_base.
	friend base_t;
	//!\endcond

public:
	/*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr cigar_op() noexcept : base_t{} {}
    constexpr cigar_op(cigar_op const &) = default;
    constexpr cigar_op(cigar_op &&) = default;
    constexpr cigar_op & operator=(cigar_op const &) = default;
    constexpr cigar_op & operator=(cigar_op &&) = default;
    ~cigar_op() = default;
    //!\}

protected:
	//!\privatesection

    //!\brief Value to char conversion table.
	static constexpr char_type rank_to_char[value_size]
	{
		'M',
		'D',
		'I',
		'S',
		'H',
		'N',
		'P',
		'X',
		'='
	};

	//!\brief Char to value conversion table.
	static constexpr std::array<rank_type, 256> char_to_rank
	{
		[] () constexpr
		{
			std::array<rank_type, 256> ret{};

			// reverse mapping for characters
			for (size_t rnk = 0u; rnk < value_size; ++rnk)
			{
				ret[rank_to_char[rnk] ] = rnk;
			}

			return ret;
		}()
	};
};

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::cigar_op char literal.
 * \relates seqan3::cigar_op
 * \returns seqan3::cigar_op
 */
cigar_op operator""_cigar_op(char const c) noexcept
{
    return cigar_op{}.assign_char_strict(c);
}
//!\}
} // namespace seqan3
