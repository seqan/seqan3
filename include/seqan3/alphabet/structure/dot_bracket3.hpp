// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides the dot bracket format for RNA structure.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/structure/concept.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dot_bracket3
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The three letter RNA structure alphabet of the characters ".()".
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
 * The brackets denote RNA base pair interactions. Every left bracket must have a corresponding right bracket.
 * Pseudoknots cannot be expressed in this format. A dot (.) represents a character that is not paired.
 *
 *```console
 *     GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
 *     (((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).
 *```
 *
 * ### Example
 *
 * \include test/snippet/alphabet/structure/dot_bracket3.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
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
    constexpr dot_bracket3() noexcept = default;                                 //!< Defaulted.
    constexpr dot_bracket3(dot_bracket3 const &) noexcept = default;             //!< Defaulted.
    constexpr dot_bracket3(dot_bracket3 &&) noexcept = default;                  //!< Defaulted.
    constexpr dot_bracket3 & operator=(dot_bracket3 const &) noexcept = default; //!< Defaulted.
    constexpr dot_bracket3 & operator=(dot_bracket3 &&) noexcept = default;      //!< Defaulted.
    ~dot_bracket3() noexcept = default;                                          //!< Defaulted.

    //!\}

    //!\name RNA structure properties
    //!\{

    /*!\brief Check whether the character represents a rightward interaction in an RNA structure.
     * \returns True if the letter represents a rightward interaction, False otherwise.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_pair_open() const noexcept
    {
        return to_rank() == 1u;
    }

    /*!\brief Check whether the character represents a leftward interaction in an RNA structure.
     * \returns True if the letter represents a leftward interaction, False otherwise.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_pair_close() const noexcept
    {
        return to_rank() == 2u;
    }

    /*!\brief Check whether the character represents an unpaired position in an RNA structure.
     * \returns True if the letter represents an unpaired site, False otherwise.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_unpaired() const noexcept
    {
        return to_rank() == 0u;
    }

    /*!\brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions, up to a certain depth.
     * \details It is the number of distinct pairs of interaction symbols the format supports. The value 1 denotes no
     * pseudoknot support.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    static constexpr uint8_t max_pseudoknot_depth{1u};

    /*!\brief Get an identifier for a pseudoknotted interaction,
     * where opening and closing brackets of the same type have the same id.
     * \returns The pseudoknot id (always 0) if alph denotes an interaction, and no value otherwise.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr std::optional<uint8_t> pseudoknot_id() const noexcept
    {
        if (is_unpaired())
            return std::nullopt;
        else
            return 0;
    }
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'.', '(', ')'};

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

    //!\copydoc seqan3::dna4::char_to_rank_table
    static constexpr std::array<rank_type, 256> char_to_rank_table{[]() constexpr
                                                                   {
                                                                       std::array<rank_type, 256> rank_table{};

                                                                       // Value-initialisation of std::array does usually initialise. `fill` is explicit.
                                                                       rank_table.fill(0u);

                                                                       // canonical
                                                                       rank_table['.'] = 0u;
                                                                       rank_table['('] = 1u;
                                                                       rank_table[')'] = 2u;

                                                                       return rank_table;
                                                                   }()};
};

inline namespace literals
{

/*!\name Structure literals
 * \{
 */
/*!\brief The seqan3::db3 char literal.
 * \relatesalso seqan3::dot_bracket3
 * \param[in] ch The character to represent as dot bracket.
 * \returns seqan3::dot_bracket3
 *
 * You can use this char literal to assign a seqan3::dot_bracket3 character:
 * \include test/snippet/alphabet/structure/dot_bracket3_char_literal.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
constexpr dot_bracket3 operator""_db3(char const ch) noexcept
{
    return dot_bracket3{}.assign_char(ch);
}

/*!\brief The seqan3::db3 string literal.
 * \relatesalso seqan3::dot_bracket3
 * \param[in] str A pointer to the character string to assign.
 * \param[in] len The size of the character string to assign.
 * \returns std::vector<seqan3::dot_bracket3>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::dot_bracket3>:
 * \include test/snippet/alphabet/structure/dot_bracket3_literal.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL std::vector<dot_bracket3> operator""_db3(char const * str, std::size_t len)
{
    std::vector<dot_bracket3> vec;
    vec.resize(len);

    for (size_t idx = 0ul; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}
//!\}

} // namespace literals

} // namespace seqan3
