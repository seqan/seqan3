// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides the WUSS format for RNA structure.
 */

#pragma once

#include <cmath>
#include <limits>
#include <vector>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/structure/concept.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

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
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \ingroup alphabet_structure
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
 * ### Example
 *
 * \include test/snippet/alphabet/structure/wuss.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
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

protected:
    using typename base_t::char_type;
    using typename base_t::rank_type;

public:
    using base_t::alphabet_size;
    using base_t::to_char;
    using base_t::to_rank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr wuss() noexcept = default;                         //!< Defaulted.
    constexpr wuss(wuss const &) noexcept = default;             //!< Defaulted.
    constexpr wuss(wuss &&) noexcept = default;                  //!< Defaulted.
    constexpr wuss & operator=(wuss const &) noexcept = default; //!< Defaulted.
    constexpr wuss & operator=(wuss &&) noexcept = default;      //!< Defaulted.
    ~wuss() noexcept = default;                                  //!< Defaulted.

    //!\}

    /*!\name RNA structure properties
     * \{
     */
    /*!\brief Check whether the character represents a rightward interaction in an RNA structure.
     * \returns True if the letter represents a rightward interaction, False otherwise.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_pair_open() const noexcept
    {
        return interaction_tab[to_rank()] < 0;
    }

    /*!\brief Check whether the character represents a leftward interaction in an RNA structure.
     * \returns True if the letter represents a leftward interaction, False otherwise.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_pair_close() const noexcept
    {
        return interaction_tab[to_rank()] > 0;
    }

    /*!\brief Check whether the character represents an unpaired position in an RNA structure.
     * \returns True if the letter represents an unpaired site, False otherwise.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_unpaired() const noexcept
    {
        return interaction_tab[to_rank()] == 0;
    }

    /*!\brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions, up to a certain depth.
     *        It is the number of distinct pairs of interaction symbols the format supports: 4..30 (depends on size)
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    // formula: (alphabet size - 7 unpaired characters) / 2, as every bracket exists as opening/closing pair
    static constexpr uint8_t max_pseudoknot_depth{static_cast<uint8_t>((alphabet_size - 7) / 2)};

    /*!\brief Get an identifier for a pseudoknotted interaction, where opening and closing brackets of the same
     *        type have the same id.
     * \returns The pseudoknot id, if alph denotes an interaction, and no value otherwise.
     * \details
     * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr std::optional<uint8_t> pseudoknot_id() const noexcept
    {
        if (interaction_tab[to_rank()] != 0)
            return std::abs(interaction_tab[to_rank()]) - 1;
        else
            return std::nullopt; // unpaired
    }
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
    }

    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr std::array<char_type, alphabet_size> rank_to_char_table{
        []() constexpr
        {
            std::array<char_type, alphabet_size>
                chars{'.', ':', ',', '-', '_', '~', ';', '<', '(', '[', '{', '>', ')', ']', '}'};

            // pseudoknot letters
            for (rank_type rnk = 15u; rnk + 1u < alphabet_size; rnk += 2u)
            {
                char_type const off = static_cast<char_type>((rnk - 15u) / 2u);
                chars[rnk] = 'A' + off;
                chars[rnk + 1u] = 'a' + off;
            }

            return chars;
        }()};

    //!\copydoc seqan3::dna4::char_to_rank_table
    static constexpr std::array<rank_type, 256> char_to_rank_table{[]() constexpr
                                                                   {
                                                                       std::array<rank_type, 256> rank_table{};

                                                                       rank_table.fill(6u);

                                                                       // set alphabet values
                                                                       for (rank_type rnk = 0u; rnk < alphabet_size;
                                                                            ++rnk)
                                                                           rank_table[rank_to_char_table[rnk]] = rnk;

                                                                       return rank_table;
                                                                   }()};

    /*!\brief Lookup table for interactions: unpaired (0), pair-open (< 0), pair-close (> 0).
     * Paired brackets have the same absolute value.
     */
    static constexpr std::array<int8_t, SIZE> interaction_tab{
        []() constexpr
        {
            static_assert(static_cast<int16_t>(std::numeric_limits<int8_t>::max()) >= SIZE);
            static_assert(-static_cast<int16_t>(std::numeric_limits<int8_t>::min()) >= SIZE);

            std::array<int8_t, alphabet_size> interaction_table{};
            int8_t cnt_open = 0;
            int8_t cnt_close = 0;

            for (rank_type rnk = 0u; rnk <= 6u; ++rnk)
                interaction_table[rnk] = 0;

            for (rank_type rnk = 7u; rnk <= 10u; ++rnk)
                interaction_table[rnk] = --cnt_open;

            for (rank_type rnk = 11u; rnk <= 14u; ++rnk)
                interaction_table[rnk] = ++cnt_close;

            for (rank_type rnk = 15u; rnk + 1u < alphabet_size; rnk += 2u)
            {
                interaction_table[rnk] = --cnt_open;
                interaction_table[rnk + 1u] = ++cnt_close;
            }

            return interaction_table;
        }()};
};

/*!\brief Alias for the default type wuss51.
 * \relates seqan3::wuss
 */
using wuss51 = wuss<51>;

inline namespace literals
{

/*!\name Structure literals
 * \{
 */
/*!\brief The seqan3::wuss51 char literal.
 * \relatesalso seqan3::wuss
 * \param[in] ch The character to represent as wuss.
 * \returns seqan3::wuss51
 *
 * You can use this char literal to assign a seqan3::wuss51 character:
 * \include test/snippet/alphabet/structure/wuss_char_literal.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
constexpr wuss51 operator""_wuss51(char const ch) noexcept
{
    return wuss51{}.assign_char(ch);
}

/*!\brief The seqan3::wuss51 string literal.
 * \relatesalso seqan3::wuss
 * \param[in] str A pointer to the character string to assign.
 * \param[in] len The size of the character string to assign.
 * \returns std::vector<seqan3::wuss51>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::wuss51>:
 * \include test/snippet/alphabet/structure/wuss_literal.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL std::vector<wuss51> operator""_wuss51(char const * str, std::size_t len)
{
    std::vector<wuss51> vec;
    vec.resize(len);

    for (size_t idx = 0ul; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}
//!\}

} // namespace literals

} // namespace seqan3
