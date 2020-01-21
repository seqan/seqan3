// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides the WUSS format for RNA structure.
 */

#pragma once

#include <cmath>
#include <vector>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/structure/concept.hpp>
#include <seqan3/core/char_operations/transform.hpp>

// ------------------------------------------------------------------
// wuss
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The WUSS structure alphabet of the characters `.<>:,-_~;()[]{}AaBbCcDd`...
 * \tparam SIZE The alphabet size defaults to 51 and must be an odd number in range 15..67.
 *              It determines the allowed pseudoknot depth by adding characters AaBb..Zz to the alphabet.
 * \implements seqan3::rna_structure_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
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
 * ###Usage
 * The following code example creates a wuss vector, modifies it, and prints the result to stderr.
 * \include test/snippet/alphabet/structure/wuss_general.cpp
 */
template <uint8_t SIZE = 51>
class wuss : public alphabet_base<wuss<SIZE>, SIZE>
{
    static_assert(SIZE >= 15 && SIZE <= 67 && SIZE % 2 == 1,
                  "The wuss<> alphabet size must be an odd number in range 15..67.");

private:
    //!\brief The base class.
    using base_t = alphabet_base<wuss<SIZE>, SIZE>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

public:
    using base_t::alphabet_size;
    using base_t::to_rank;
    using base_t::to_char;
    using typename base_t::rank_type;
    using typename base_t::char_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr wuss()                         noexcept = default; //!< Defaulted.
    constexpr wuss(wuss const &)             noexcept = default; //!< Defaulted.
    constexpr wuss(wuss &&)                  noexcept = default; //!< Defaulted.
    constexpr wuss & operator=(wuss const &) noexcept = default; //!< Defaulted.
    constexpr wuss & operator=(wuss &&)      noexcept = default; //!< Defaulted.
    ~wuss()                                  noexcept = default; //!< Defaulted.
    //!\}

    /*!\name RNA structure properties
     * \{
     */
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
    static constexpr uint8_t max_pseudoknot_depth{static_cast<uint8_t>((alphabet_size - 7) / 2)};

    /*!\brief Get an identifier for a pseudoknotted interaction,
     * where opening and closing brackets of the same type have the same id.
     * \returns The pseudoknot id, if alph denotes an interaction, and no value otherwise.
     * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
     */
    constexpr std::optional<uint8_t> pseudoknot_id() const noexcept
    {
        if (interaction_tab[to_rank()] != 0)
            return std::abs(interaction_tab[to_rank()]) - 1;
        else
            return std::nullopt; // unpaired
    }
    //!\}

protected:
    //!\privatesection
    //!\brief Value-to-char conversion table.
    static constexpr std::array<char_type, alphabet_size> rank_to_char
    {
        [] () constexpr
        {
            std::array<char_type, alphabet_size> chars
            {
                '.', ':', ',', '-', '_', '~', ';', '<', '(', '[', '{', '>', ')', ']', '}'
            };

            // pseudoknot letters
            for (rank_type rnk = 15u; rnk + 1u < alphabet_size; rnk += 2u)
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
                rnk = 6u;

            // set alphabet values
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
                rank_table[rank_to_char[rnk]] = rnk;
            return rank_table;
        } ()
    };

    /*!\brief Lookup table for interactions: unpaired (0), pair-open (< 0), pair-close (> 0).
     * Paired brackets have the same absolute value.
     */
    static std::array<int8_t, SIZE> const interaction_tab;
};

template <uint8_t SIZE>
constexpr std::array<int8_t, SIZE> wuss<SIZE>::interaction_tab = [] () constexpr
{
    std::array<int8_t, alphabet_size> interaction_table{};
    int cnt_open = 0;
    int cnt_close = 0;

    for (rank_type rnk = 0u; rnk <= 6u; ++rnk)
    {
        interaction_table[rnk] = 0;
    }

    for (rank_type rnk = 7u; rnk <= 10u; ++rnk)
    {
        interaction_table[rnk] = --cnt_open;
    }

    for (rank_type rnk = 11u; rnk <= 14u; ++rnk)
    {
        interaction_table[rnk] = ++cnt_close;
    }

    for (rank_type rnk = 15u; rnk + 1u < alphabet_size; rnk += 2u)
    {
        interaction_table[rnk]      = --cnt_open;
        interaction_table[rnk + 1u] = ++cnt_close;
    }

    return interaction_table;
} ();

//!\brief Alias for the default type wuss51.
//!\relates seqan3::wuss
using wuss51 = wuss<51>;

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::wuss51 string literal.
 * \relates seqan3::wuss
 * \param[in] str A pointer to the character string to assign.
 * \param[in] len The size of the character string to assign.
 * \returns std::vector<seqan3::wuss51>
 *
 * You can use this string literal to easily assign to a vector of seqan3::wuss51 characters:
 * \include test/snippet/alphabet/structure/wuss_literal.cpp
 */
inline std::vector<wuss51> operator""_wuss51(const char * str, std::size_t len)
{
    std::vector<wuss51> vec;
    vec.resize(len);

    for (size_t idx = 0ul; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}

/*!\brief The seqan3::wuss51 char literal.
 * \relates seqan3::wuss
 * \param[in] ch The character to represent as wuss.
 * \returns seqan3::wuss51
 *
 * You can use this string literal to assign a seqan3::wuss51 character.
 * For different wuss alphabet sizes the `assign_char` function must be used.
 * \include test/snippet/alphabet/structure/wuss_char_literal.cpp
 */
constexpr wuss51 operator""_wuss51(char const ch) noexcept
{
    return wuss51{}.assign_char(ch);
}

//!\}

} // namespace seqan3
