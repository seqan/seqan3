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
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains seqan3::aa27, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/aminoacid/aminoacid_base.hpp>
#include <seqan3/io/stream/char_operations.hpp>

namespace seqan3
{
/*!\brief The twenty-seven letter amino acid alphabet.
 * \ingroup aminoacid
 * \implements seqan3::aminoacid_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 * \details
 * The alphabet consists of letters A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X,
 * Y, Z, *
 *
 * The alphabet may be brace initialized from the static letter members. Note that you cannot
 * assign the alphabet by using letters of type `char`, but you instead have to use the
 * function seqan3::aa27::assign_char().
 *
 * \snippet test/snippet/alphabet/aminoacid/aa27.cpp construction
 */

class aa27 : public aminoacid_base<aa27, 27>
{
private:
    //!\brief The base class.
    using base_t = aminoacid_base<aa27, 27>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aa27() : base_t{} {}
    constexpr aa27(aa27 const &) = default;
    constexpr aa27(aa27 &&) = default;
    constexpr aa27 & operator=(aa27 const &) = default;
    constexpr aa27 & operator=(aa27 &&) = default;
    ~aa27() = default;

    using base_t::base_t;
    //!\}

protected:
    //!\brief Value to char conversion table.
    static constexpr char_type rank_to_char[value_size]
    {
        'A',
        'B',
        'C',
        'D',
        'E',
        'F',
        'G',
        'H',
        'I',
        'J',
        'K',
        'L',
        'M',
        'N',
        'O',
        'P',
        'Q',
        'R',
        'S',
        'T',
        'U',
        'V',
        'W',
        'X',
        'Y',
        'Z',
        '*'
    };

    //!\brief Char to value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = 23; // value of 'X'

            // reverse mapping for characters and their lowercase
            for (rank_type rnk = 0u; rnk < value_size; ++rnk)
            {
                ret[static_cast<rank_type>(         rank_to_char[rnk]) ] = rnk;
                ret[static_cast<rank_type>(to_lower(rank_to_char[rnk]))] = rnk;
            }

            return ret;
        }()
    };
};

} // namespace seqan3

// ------------------------------------------------------------------
// metafunctions
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Helper metafunction that identifies aa27 as an amino acid alphabet.
//!\ingroup aminoacid
template <>
struct is_aminoacid<aa27> : std::true_type {};

} // namespace seqan3

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{
//!\brief Alias for an std::vector of seqan3::aa27.
//!\relates aa27
using aa27_vector = std::vector<aa27>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3
{

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::aa27 char literal.
 * \param[in] c The character to assign.
 * \relates seqan3::aa27
 * \returns seqan3::aa27
 *
 * \snippet test/snippet/alphabet/aminoacid/aa27.cpp char_literal
 *
 */
constexpr aa27 operator""_aa27(char const c) noexcept
{
    return aa27{}.assign_char(c);
}

/*!\brief The seqan3::aa27 string literal.
 * \param[in] s A pointer to the character string to assign.
 * \param[in] n The size of the character string to assign.
 * \relates seqan3::aa27
 * \returns seqan3::aa27_vector
 *
 * You can use this string literal to easily assign to aa27_vector:
 *
 * \snippet test/snippet/alphabet/aminoacid/aa27.cpp literal
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3!
 */

inline aa27_vector operator""_aa27(const char * s, std::size_t n)
{
    aa27_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace seqan3
