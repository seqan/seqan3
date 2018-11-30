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

#include <cmath>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/io/stream/char_operations.hpp>

// ------------------------------------------------------------------
// wuss
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The WUSS structure alphabet of the characters `.<>:,-_~;()[]{}AaBbCcDd`...
 * \tparam SIZE The alphabet size defaults to 50 and must be an odd number in range 15..67.
 *              It determines the allowed pseudoknot depth by adding characters AaBb..Zz to the alphabet.
 * \implements seqan3::rna_structure_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
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
 * \snippet test/snippet/alphabet/structure/wuss.cpp general
 */
template <uint8_t SIZE = 51>
class wuss : public alphabet_base<wuss<SIZE>, SIZE>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<wuss<SIZE>, SIZE>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

    //!\brief Required for deferred instantiation of member "enum" values.
    using this_type_deferred = typename base_t::derived_t;

public:
    using base_t::value_size;
    using base_t::to_rank;
    using base_t::to_char;
    using typename base_t::rank_type;
    using typename base_t::char_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr wuss() : base_t{} {}
    constexpr wuss(wuss const &) = default;
    constexpr wuss(wuss &&) = default;
    constexpr wuss & operator=(wuss const &) = default;
    constexpr wuss & operator=(wuss &&) = default;
    ~wuss() = default;
    //!\}

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface. *Don't worry about the `internal_type`.*
     * The pseudoknot letters are not accessible in this way, please use the string literals.
     */
    //!\{

    //!\brief `.` not paired (insertion to known structure)
    static this_type_deferred constexpr UNPAIRED  = this_type_deferred{}.assign_char('.');
    //!\brief `:` not paired (external residue outside structure)
    static this_type_deferred constexpr UNPAIRED1 = this_type_deferred{}.assign_char(':');
    //!\brief `,` not paired (multifurcation loop)
    static this_type_deferred constexpr UNPAIRED2 = this_type_deferred{}.assign_char(',');
    //!\brief `-` not paired (bulge, interior loop)
    static this_type_deferred constexpr UNPAIRED3 = this_type_deferred{}.assign_char('-');
    //!\brief `_` not paired (hairpin loop)
    static this_type_deferred constexpr UNPAIRED4 = this_type_deferred{}.assign_char('_');
    //!\brief `~` not paired (due to local alignment)
    static this_type_deferred constexpr UNPAIRED5 = this_type_deferred{}.assign_char('~');
    //!\brief `;` not paired
    static this_type_deferred constexpr UNPAIRED6 = this_type_deferred{}.assign_char(';');

    //!\brief `<` bracket left (simple terminal stem)
    static this_type_deferred constexpr PAIR_OPEN   = this_type_deferred{}.assign_char('<');
    //!\brief `(` bracket left (internal helix enclosing `<>`)
    static this_type_deferred constexpr PAIR_OPEN1  = this_type_deferred{}.assign_char('(');
    //!\brief `[` bracket left (internal helix enclosing `()`)
    static this_type_deferred constexpr PAIR_OPEN2  = this_type_deferred{}.assign_char('[');
    //!\brief `{` bracket left (internal helix enclosing `[]`)
    static this_type_deferred constexpr PAIR_OPEN3  = this_type_deferred{}.assign_char('{');

    //!\brief `>` bracket right (simple terminal stem)
    static this_type_deferred constexpr PAIR_CLOSE  = this_type_deferred{}.assign_char('>');
    //!\brief `)` bracket right (internal helix enclosing `<>`)
    static this_type_deferred constexpr PAIR_CLOSE1 = this_type_deferred{}.assign_char(')');
    //!\brief `]` bracket right (internal helix enclosing `()`)
    static this_type_deferred constexpr PAIR_CLOSE2 = this_type_deferred{}.assign_char(']');
    //!\brief `}` bracket right (internal helix enclosing `[]`)
    static this_type_deferred constexpr PAIR_CLOSE3 = this_type_deferred{}.assign_char('}');
    // pseudoknots not accessible
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
    static constexpr uint8_t max_pseudoknot_depth{(value_size - 7) / 2};

    /*!\brief Get an identifier for a pseudoknotted interaction.
     * Opening and closing brackets of the same type have the same id.
     * \returns The pseudoknot id, if alph denotes an interaction, and no value otherwise.
     * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
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
    //!\brief Value-to-char conversion table.
    static constexpr std::array<char_type, value_size> rank_to_char
    {
        [] () constexpr
        {
            std::array<char_type, value_size> chars
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
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> rank_table{};

            // initialize with unpaired (std::array::fill unfortunately not constexpr)
            for (rank_type & rnk : rank_table)
                rnk = 6; // ::UNPAIRED6;

            // set alphabet values
            for (rank_type rnk = 0u; rnk < value_size; ++rnk)
                rank_table[rank_to_char[rnk]] = rnk;
            return rank_table;
        } ()
    };

    //!\brief Lookup table for interactions: unpaired (1), pair-open (2), pair-close (3).
    static std::array<int8_t, SIZE> const interaction_tab;
};

template <uint8_t SIZE>
constexpr std::array<int8_t, SIZE> wuss<SIZE>::interaction_tab = [] () constexpr
{
    std::array<int8_t, value_size> interaction_table{};
    int cnt_open = 0;
    int cnt_close = 0;

    for (rank_type rnk = UNPAIRED.to_rank();
            rnk <= UNPAIRED6.to_rank();
            ++rnk)
    {
        interaction_table[rnk] = 0;
    }

    for (rank_type rnk = PAIR_OPEN.to_rank();
            rnk <= PAIR_OPEN3.to_rank();
            ++rnk)
    {
        interaction_table[rnk] = --cnt_open;
    }

    for (rank_type rnk = PAIR_CLOSE.to_rank();
            rnk <= PAIR_CLOSE3.to_rank();
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
} ();

//!\brief Alias for the default type wuss51.
typedef wuss<51> wuss51;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief wuss literal
 * \relates seqan3::wuss
 * \returns std::vector<seqan3::wuss51>
 *
 * You can use this string literal to easily assign to a vector of wuss characters:
 *
 *```.cpp
 *     std::vector<wuss<>> foo{".<..>."_wuss51};
 *     std::vector<wuss<>> bar = ".<..>."_wuss51;
 *     auto bax = ".<..>."_wuss51;
 *```
 */
inline std::vector<wuss51> operator""_wuss51(const char * str, std::size_t len)
{
    std::vector<wuss51> vec;
    vec.resize(len);

    for (size_t idx = 0; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}

} // namespace seqan3
