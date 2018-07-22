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
 * \brief Contains the WUSS format for RNA structure.
 */

#pragma once

#include <cassert>
#include <cmath>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

// ------------------------------------------------------------------
// wuss
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The WUSS structure alphabet of the characters `.<>:,-_~;()[]{}AaBbCcDd`...
 * \tparam SIZE The alphabet size defaults to 50 and must be an odd number in range 15..67.
 *              It determines the allowed pseudoknot depth by adding characters AaBb..Zz to the alphabet.
 * \implements seqan3::rna_structure_concept
 * \ingroup structure
 *
 * \details
 * The symbols `.:,-_~;` denote unpaired characters, brackets `<>()[]{}` represent base pair interactions and
 * `AaBbCcDd`... form pseudoknots in the structure. The default alphabet has size 51 (letters until `Rr`).
 * The size can be varied with the optional template parameter between 15 (no letters for pseudoknots) and 67
 * (all `Aa`-`Zz` for pseudoknots).
 *
 *```console
 * <<<___>>>,,<<<__>>>
 * <<<<_AAAA____>>>>aaaa
 *```
 *
 * \par Usage
 * The following code example creates a wuss vector, modifies it, and prints the result to stdout.
 * ```cpp
 *     // create vector
 *     std::vector<wuss51> vec{wuss51::UNPAIRED, wuss51::PAIR_CLOSE, wuss51::PAIR_CLOSE};
 *     // modify and print
 *     vec[1] = wuss51::PAIR_OPEN;
 *     for (wuss51 chr : vec)
 *         std::cout << chr;  // .<>
 * ```
 */
template <uint8_t SIZE = 51>
struct wuss
{
    // check alphabet size constraint
    static_assert(SIZE >= 15 && SIZE <= 67 && SIZE % 2 == 1);

    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = uint8_t;

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface. *Don't worry about the `internal_type`.*
     * The pseudoknot letters are not accessible in this way, please use the string literals.
     */
    //!\{

    //!\brief `.` not paired (insertion to known structure)
    static const wuss UNPAIRED;
    //!\brief `:` not paired (external residue outside structure)
    static const wuss UNPAIRED1;
    //!\brief `,` not paired (multifurcation loop)
    static const wuss UNPAIRED2;
    //!\brief `-` not paired (bulge, interior loop)
    static const wuss UNPAIRED3;
    //!\brief `_` not paired (hairpin loop)
    static const wuss UNPAIRED4;
    //!\brief `~` not paired (due to local alignment)
    static const wuss UNPAIRED5;
    //!\brief `;` not paired
    static const wuss UNPAIRED6;

    //!\brief `<` bracket left (simple terminal stem)
    static const wuss PAIR_OPEN;
    //!\brief `(` bracket left (internal helix enclosing `<>`)
    static const wuss PAIR_OPEN1;
    //!\brief `[` bracket left (internal helix enclosing `()`)
    static const wuss PAIR_OPEN2;
    //!\brief `{` bracket left (internal helix enclosing `[]`)
    static const wuss PAIR_OPEN3;

    //!\brief `>` bracket right (simple terminal stem)
    static const wuss PAIR_CLOSE;
    //!\brief `)` bracket right (internal helix enclosing `<>`)
    static const wuss PAIR_CLOSE1;
    //!\brief `]` bracket right (internal helix enclosing `()`)
    static const wuss PAIR_CLOSE2;
    //!\brief `}` bracket right (internal helix enclosing `[]`)
    static const wuss PAIR_CLOSE3;
    // pseudoknots not accessible
    //!\}

    //!\name Read functions
    //!\{

    /*!
     * \brief Get the letter as a character of char_type.
     * \returns The character representation of this wuss letter.
     */
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    /*!\brief Get the letter's numeric value or rank in the alphabet.
     * \returns The numeric representation of this wuss letter.
     */
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }
    //!\}

    //!\name Write functions
    //!\{

    /*!
     * \brief Assign from a character.
     * \param chr The character that is assigned.
     * \returns The resulting wuss character.
     */
    constexpr wuss & assign_char(char_type const chr) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        _value = char_to_value[static_cast<index_t>(chr)];
        return *this;
    }

    /*!\brief Assign from a numeric value.
     * \param rnk The rank value that is assigned.
     * \returns The resulting wuss character.
     */
    constexpr wuss & assign_rank(rank_type const rnk)
    {
        assert(rnk < value_size);
        _value = static_cast<internal_type>(rnk);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{SIZE};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(wuss const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(wuss const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(wuss const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(wuss const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(wuss const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(wuss const & rhs) const noexcept
    {
        return _value >= rhs._value;
    }
    //!\}

    //!\name RNA structure properties
    //!\{

    /*!\brief Check whether the character represents a rightward interaction in an RNA structure.
     * \returns True if the letter represents a rightward interaction, False otherwise.
     */
    constexpr bool is_pair_open() const noexcept
    {
        return interaction_tab[to_rank()] < 0;
    }

    /*!\brief Check whether the character represents a leftward interaction in an RNA structure.
     * \returns True if the letter represents a leftward interaction, False otherwise.
     */
    constexpr bool is_pair_close() const noexcept
    {
        return interaction_tab[to_rank()] > 0;
    }

    /*!\brief Check whether the character represents an unpaired position in an RNA structure.
     * \returns True if the letter represents an unpaired site, False otherwise.
     */
    constexpr bool is_unpaired() const noexcept
    {
        return interaction_tab[to_rank()] == 0;
    }

    /*!\brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions, up to a certain depth.
     *        It is the number of distinct pairs of interaction symbols the format supports: 4..30 (depends on size)
     */
    // formula: (alphabet size - 7 unpaired characters) / 2, as every bracket exists as opening/closing pair
    static constexpr uint8_t pseudoknot_support{(value_size - 7) / 2};

    /*!\brief Get an identifier for a pseudoknotted interaction.
     * Opening and closing brackets of the same type have the same id.
     * \returns The pseudoknot id, if alph denotes an interaction, and no value otherwise.
     * It is guaranteed to be smaller than seqan3::pseudoknot_support.
     */
    constexpr std::optional<uint8_t> pseudoknot_id() const noexcept
    {
        if (interaction_tab[to_rank()] != 0)
            return std::abs(interaction_tab[to_rank()]) - 1;
        else
            return std::nullopt;
    }
    //!\}

protected:
    //!\privatesection
    /*!\brief The internal type is a strictly typed enum.
     *
     * This is done to prevent aggregate initialization from numbers and/or chars.
     * It is has the drawback that it also introduces a scope which in turn makes
     * the static "letter values " members necessary.
     */
    enum struct internal_type : rank_type
    {
        UNPAIRED,    // not paired .
        UNPAIRED1,   // not paired :
        UNPAIRED2,   // not paired ,
        UNPAIRED3,   // not paired -
        UNPAIRED4,   // not paired _
        UNPAIRED5,   // not paired ~
        UNPAIRED6,   // not paired ;
        PAIR_OPEN,   // bracket left <
        PAIR_OPEN1,  // bracket left (
        PAIR_OPEN2,  // bracket left [
        PAIR_OPEN3,  // bracket left {
        PAIR_CLOSE,  // bracket right >
        PAIR_CLOSE1, // bracket right )
        PAIR_CLOSE2, // bracket right ]
        PAIR_CLOSE3  // bracket right }
    };

    //!\brief Value-to-char conversion table.
    static constexpr std::array<char_type, value_size> value_to_char
    {
        [] () constexpr
        {
            std::array<char_type , value_size> chars
            {
                '.', ':', ',', '-', '_', '~', ';', '<', '(', '[', '{', '>', ')', ']', '}'
            };

            // pseudoknot letters
            for (rank_type rnk = 15u; rnk + 1u < value_size; rnk += 2u)
            {
                char_type const off = static_cast<char_type>((rnk - 15u) / 2u);
                chars[rnk] = 'A' + off;
                chars[rnk + 1u] = 'a' + off;
            }

            return chars;
        } ()
    };

    //!\brief Char-to-value conversion table.
    static constexpr std::array<internal_type, 256> char_to_value
    {
        [] () constexpr
        {
            std::array<internal_type, 256> rank_table{};

            // initialize with unpaired (std::array::fill unfortunately not constexpr)
            for (internal_type & rnk : rank_table)
                rnk = internal_type::UNPAIRED6;

            // set alphabet values
            for (rank_type rnk = 0u; rnk < value_size; ++rnk)
                rank_table[value_to_char[rnk]] = static_cast<internal_type>(rnk);
            return rank_table;
        } ()
    };

    //!\brief Lookup table for interactions: unpaired (1), pair-open (2), pair-close (3).
    static constexpr std::array<int8_t, value_size> interaction_tab
    {
        [] () constexpr
        {
            std::array<int8_t, value_size> interaction_table{};
            int cnt_open = 0;
            int cnt_close = 0;

            for (rank_type rnk = static_cast<rank_type>(internal_type::UNPAIRED);
                 rnk <= static_cast<rank_type>(internal_type::UNPAIRED6);
                 ++rnk)
            {
                interaction_table[rnk] = 0;
            }

            for (rank_type rnk = static_cast<rank_type>(internal_type::PAIR_OPEN);
                 rnk <= static_cast<rank_type>(internal_type::PAIR_OPEN3);
                 ++rnk)
            {
                interaction_table[rnk] = --cnt_open;
            }

            for (rank_type rnk = static_cast<rank_type>(internal_type::PAIR_CLOSE);
                 rnk <= static_cast<rank_type>(internal_type::PAIR_CLOSE3);
                 ++rnk)
            {
                interaction_table[rnk] = ++cnt_close;
            }

            for (rank_type rnk = 15u; rnk + 1 < value_size; rnk += 2u)
            {
                interaction_table[rnk]     = --cnt_open;
                interaction_table[rnk + 1] = ++cnt_close;
            }

            return interaction_table;
        } ()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::UNPAIRED{internal_type::UNPAIRED};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::UNPAIRED1{internal_type::UNPAIRED1};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::UNPAIRED2{internal_type::UNPAIRED2};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::UNPAIRED3{internal_type::UNPAIRED3};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::UNPAIRED4{internal_type::UNPAIRED4};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::UNPAIRED5{internal_type::UNPAIRED5};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::UNPAIRED6{internal_type::UNPAIRED6};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_OPEN{internal_type::PAIR_OPEN};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_OPEN1{internal_type::PAIR_OPEN1};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_OPEN2{internal_type::PAIR_OPEN2};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_OPEN3{internal_type::PAIR_OPEN3};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_CLOSE{internal_type::PAIR_CLOSE};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_CLOSE1{internal_type::PAIR_CLOSE1};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_CLOSE2{internal_type::PAIR_CLOSE2};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::PAIR_CLOSE3{internal_type::PAIR_CLOSE3};

//!\brief Alias for the default type wuss51.
typedef wuss<51> wuss51;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief wuss literal
 * \relates seqan3::wuss
 * \returns std::vector<seqan3::wuss51>
 *
 * You can use this string literal to easily assign to a vector of wuss characters:
 *
 *```.cpp
 *     using namespace seqan3::literal;
 *     std::vector<wuss<>> foo{".<..>."_wuss51};
 *     std::vector<wuss<>> bar = ".<..>."_wuss51;
 *     auto bax = ".<..>."_wuss51;
 *```
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */
inline std::vector<wuss51> operator""_wuss51(const char * str, std::size_t len)
{
    std::vector<wuss51> vec;
    vec.resize(len);

    for (size_t idx = 0; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}

} // namespace seqan3::literal
