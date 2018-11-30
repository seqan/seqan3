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
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Contains the dssp format for protein structure.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/io/stream/char_operations.hpp>

// ------------------------------------------------------------------
// dssp9
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The protein structure alphabet of the characters "HGIEBTSCX".
 * \implements seqan3::alphabet_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 * \ingroup structure
 *
 * \details
 * The DSSP annotation links structure elements to protein sequences.
 * Originally created with 7 letters as a file format for the DSSP program (http://www.cmbi.ru.nl/dssp.html),
 * it is also used in the stockholm file format for structure alignments, extended by the characters C and X
 * (https://en.wikipedia.org/wiki/Stockholm_format).
 *
 * The letter abbreviations are as follows:
 *
 * H = alpha helix
 * B = beta bridge
 * E = strand
 * G = helix-3
 * I = helix-5
 * T = turn
 * S = bend
 * C = coil/loop
 * X = unknown
 *
 * \par Usage
 * The following code example creates a dssp9 vector, modifies it, and prints the result to stdout.
 * \snippet test/snippet/alphabet/structure/dssp9.cpp general
 */
class dssp9 : public alphabet_base<dssp9, 9>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<dssp9, 9>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dssp9() : base_t{} {}
    constexpr dssp9(dssp9 const &) = default;
    constexpr dssp9(dssp9 &&) = default;
    constexpr dssp9 & operator=(dssp9 const &) = default;
    constexpr dssp9 & operator=(dssp9 &&) = default;
    ~dssp9() = default;
    //!\}

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface. *Don't worry about the `internal_type`.*
     */
    //!\{
    static const dssp9 H;
    static const dssp9 B;
    static const dssp9 E;
    static const dssp9 G;
    static const dssp9 I;
    static const dssp9 T;
    static const dssp9 S;
    static const dssp9 C;
    static const dssp9 X;
    //!\}

protected:
    //!\privatesection

    //!\brief Value-to-char conversion table.
    static constexpr char_type rank_to_char[value_size]
    {
        'H', 'B', 'E', 'G', 'I', 'T', 'S', 'C', 'X'
    };

    //!\brief Char-to-value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with X (std::array::fill unfortunately not constexpr)
            for (rank_type & rnk : ret)
                rnk = 8;

            // reverse mapping for characters
            for (rank_type rnk = 0u; rnk < value_size; ++rnk)
            {
                ret[static_cast<rank_type>(rank_to_char[rnk])] = rnk;
            }

            return ret;
        } ()
    };
};

constexpr dssp9 dssp9::H = dssp9{}.assign_char('H');
constexpr dssp9 dssp9::B = dssp9{}.assign_char('B');
constexpr dssp9 dssp9::E = dssp9{}.assign_char('E');
constexpr dssp9 dssp9::G = dssp9{}.assign_char('G');
constexpr dssp9 dssp9::I = dssp9{}.assign_char('I');
constexpr dssp9 dssp9::T = dssp9{}.assign_char('T');
constexpr dssp9 dssp9::S = dssp9{}.assign_char('S');
constexpr dssp9 dssp9::C = dssp9{}.assign_char('C');
constexpr dssp9 dssp9::X = dssp9{}.assign_char('X');

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief dssp9 literal
 * \relates seqan3::dssp9
 * \returns std::vector<seqan3::dssp9>
 *
 * You can use this string literal to easily assign to a vector of dssp9 characters:
 *
 *```.cpp
 *     std::vector<dssp9> foo{"EHHHHT"_dssp9};
 *     std::vector<dssp9> bar = "EHHHHT"_dssp9;
 *     auto bax = "EHHHHT"_dssp9;
 *```
 */
inline std::vector<dssp9> operator""_dssp9(const char * str, std::size_t len)
{
    std::vector<dssp9> vec;
    vec.resize(len);

    for (size_t idx = 0u; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}

} // namespace seqan3
