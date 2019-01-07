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
 * \brief Contains the dot bracket format for RNA structure.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/alphabet/structure/rna_structure_concept.hpp>
#include <seqan3/io/stream/char_operations.hpp>

// ------------------------------------------------------------------
// dot_bracket3
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The three letter RNA structure alphabet of the characters ".()".
 * \implements seqan3::rna_structure_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 * \ingroup structure
 *
 * \details
 * The brackets denote RNA base pair interactions. Every left bracket must have a corresponding right bracket.
 * Pseudoknots cannot be expressed in this format. A dot (.) represents a character that is not paired.
 *
 *```console
 *     GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
 *     (((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).
 *```
 *
 * \par Usage
 * The following code example creates a dot_bracket3 vector, modifies it, and prints the result to stderr.
 * \snippet test/snippet/alphabet/structure/dot_bracket3.cpp general
 */
class dot_bracket3 : public alphabet_base<dot_bracket3, 3>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<dot_bracket3, 3>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dot_bracket3() : base_t{} {}
    constexpr dot_bracket3(dot_bracket3 const &) = default;
    constexpr dot_bracket3(dot_bracket3 &&) = default;
    constexpr dot_bracket3 & operator=(dot_bracket3 const &) = default;
    constexpr dot_bracket3 & operator=(dot_bracket3 &&) = default;
    ~dot_bracket3() = default;
    //!\}

    //!\name RNA structure properties
    //!\{

    /*!\brief Check whether the character represents a rightward interaction in an RNA structure.
     * \returns True if the letter represents a rightward interaction, False otherwise.
     */
    constexpr bool is_pair_open() const noexcept
    {
        return this->to_char() == '(';
    }

    /*!\brief Check whether the character represents a leftward interaction in an RNA structure.
     * \returns True if the letter represents a leftward interaction, False otherwise.
     */
    constexpr bool is_pair_close() const noexcept
    {
        return this->to_char() == ')';
    }

    /*!\brief Check whether the character represents an unpaired position in an RNA structure.
     * \returns True if the letter represents an unpaired site, False otherwise.
     */
    constexpr bool is_unpaired() const noexcept
    {
        return this->to_char() == '.';
    }

    /*!\brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions, up to a certain depth.
     * \details It is the number of distinct pairs of interaction symbols the format supports. The value 1 denotes no
     * pseudoknot support.
     */
    static constexpr uint8_t max_pseudoknot_depth{1};
    //!\}

protected:
    //!\privatesection

    //!\brief Value-to-char conversion table.
    static constexpr char_type rank_to_char[value_size]
    {
        '.',
        '(',
        ')'
    };

    //!\brief Char-to-value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> rank_table{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (rank_type & rnk : rank_table)
                rnk = 0;

            // canonical
            rank_table['.'] = 0;
            rank_table['('] = 1;
            rank_table[')'] = 2;

            return rank_table;
        } ()
    };
};

/*!\name Literals
 * \{
 *
 * \brief The seqan3::db3 string literal.
 * \relates seqan3::dot_bracket3
 * \param[in] str A pointer to the character string to assign.
 * \param[in] len The size of the character string to assign.
 * \returns std::vector<seqan3::dot_bracket3>
 *
 * You can use this string literal to easily assign to a vector of seqan3::dot_bracket3 characters:
 * \snippet test/snippet/alphabet/structure/dot_bracket3.cpp string_literal
 */
inline std::vector<dot_bracket3> operator""_db3(const char * str, std::size_t len)
{
    std::vector<dot_bracket3> vec;
    vec.resize(len);

    for (size_t idx = 0u; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}

/*!
 * \brief The seqan3::db3 char literal.
 * \relates seqan3::dot_bracket3
 * \param[in] ch The character to represent as dot bracket.
 * \returns seqan3::dot_bracket3
 *
 * You can use this string literal to assign a seqan3::dot_bracket3 character:
 * \snippet test/snippet/alphabet/structure/dot_bracket3.cpp char_literal
 */
constexpr dot_bracket3 operator""_db3(char const ch) noexcept
{
    return dot_bracket3{}.assign_char(ch);
}

//!\}

} // namespace seqan3
