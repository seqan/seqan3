// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::dna4, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>

// ------------------------------------------------------------------
// dna4
// ------------------------------------------------------------------

namespace seqan3
{

class rna4;

/*!\brief The four letter DNA alphabet of A,C,G,T.
 * \ingroup alphabet_nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * Note that you can assign 'U' as a character to dna4 and it will silently
 * be converted to 'T'.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal \ref seqan3_dna4_char_literal "'A'_dna4" or use the
 * function seqan3::dna4::assign_char().
 *
 * \include test/snippet/alphabet/nucleotide/dna4.cpp
 *
 * If the special char conversion of IUPAC characters is not your desired behaviour, refer to our cookbook for an
 * example of \ref cookbook_custom_dna4_alphabet to change the conversion behaviour.
 *
 * \stableapi{Since version 3.1.}
 */
class dna4 : public nucleotide_base<dna4, 4>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna4, 4>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond
    //!\brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond
    //!\brief Befriend seqan3::rna4 so it can copy #char_to_rank.
    friend rna4;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna4() noexcept = default;                         //!< Defaulted.
    constexpr dna4(dna4 const &) noexcept = default;             //!< Defaulted.
    constexpr dna4(dna4 &&) noexcept = default;                  //!< Defaulted.
    constexpr dna4 & operator=(dna4 const &) noexcept = default; //!< Defaulted.
    constexpr dna4 & operator=(dna4 &&) noexcept = default;      //!< Defaulted.
    ~dna4() noexcept = default;                                  //!< Defaulted.

    using base_t::base_t;

    /*!\brief Allow implicit construction from seqan3::rna4 of the same size.
     * \details
     *
     * \include{doc} doc/fragments/dna4_implicit_conversion_from_rna4.md
     *
     * \stableapi{Since version 3.1.}
     */
    template <std::same_as<rna4> t> // template parameter t to accept incomplete type
    constexpr dna4(t const & r) noexcept
    {
        assign_rank(r.to_rank());
    }
    //!\}

    //!\brief Returns the complement of the current nucleotide.
    constexpr dna4 complement() const noexcept
    {
        return dna4{}.assign_rank(to_rank() ^ 0b11);
    }

private:
    /*!\brief The lookup table used in #rank_to_char.
     * \details
     * We would have defined these lookup tables directly within their respective constexpr functions, but at the time
     * of writing this, gcc did not (clang >= 4 did!) auto-generate lookup tables.
     *
     * ```cpp
     * static constexpr char_type rank_to_char(rank_type const rank)
     * {
     *     // not possible because of static not being allowed within a constexpr function
     *     static constexpr lookup_table = ...;
     *     return lookup_table[rank];
     * }
     *
     * static constexpr char_type rank_to_char(rank_type const rank)
     * {
     *     // up-to the compiler to optimise, no guarantee that a lookup table is used.
     *     constexpr lookup_table = ...;
     *     return lookup_table[rank];
     * }
     * ```
     *
     * \sa https://gcc.gnu.org/bugzilla/show_bug.cgi?id=99320 for the progress on gcc
     */
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'C', 'G', 'T'};

    //!\brief The rank complement table.
    static constexpr rank_type rank_complement_table[alphabet_size]{
        3, // T is complement of 'A'_dna4
        2, // G is complement of 'C'_dna4
        1, // C is complement of 'G'_dna4
        0  // A is complement of 'T'_dna4
    };

    /*!\brief Returns the complement by rank.
     * \details
     * This function is required by seqan3::nucleotide_base.
     */
    static constexpr rank_type rank_complement(rank_type const rank)
    {
        return rank_complement_table[rank];
    }

    /*!\brief Returns the character representation of rank.
     * \details
     * This function is required by seqan3::alphabet_base.
     */
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    /*!\brief Returns the rank representation of character.
     * \details
     * This function is required by seqan3::alphabet_base.
     */
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
    }

    /*!\brief The lookup table used in #char_to_rank.
     * \copydetails seqan3::dna4::rank_to_char_table
     */
    static constexpr std::array<rank_type, 256> char_to_rank_table{[]() constexpr
                                                                   {
                                                                       std::array<rank_type, 256> ret{};

                                                                       // reverse mapping for characters and their lowercase
                                                                       for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
                                                                       {
                                                                           ret[rank_to_char_table[rnk]] = rnk;
                                                                           ret[to_lower(rank_to_char_table[rnk])] = rnk;
                                                                       }

                                                                       // set U equal to T
                                                                       ret['U'] = ret['T'];
                                                                       ret['u'] = ret['t'];

                                                                       // iupac characters get special treatment, because there is no N
                                                                       ret['R'] = ret['A'];
                                                                       ret['r'] = ret['A']; // A or G
                                                                       ret['Y'] = ret['C'];
                                                                       ret['y'] = ret['C']; // C or T
                                                                       ret['S'] = ret['C'];
                                                                       ret['s'] = ret['C']; // C or G
                                                                       ret['W'] = ret['A'];
                                                                       ret['w'] = ret['A']; // A or T
                                                                       ret['K'] = ret['G'];
                                                                       ret['k'] = ret['G']; // G or T
                                                                       ret['M'] = ret['A'];
                                                                       ret['m'] = ret['A']; // A or T
                                                                       ret['B'] = ret['C'];
                                                                       ret['b'] = ret['C']; // C or G or T
                                                                       ret['D'] = ret['A'];
                                                                       ret['d'] = ret['A']; // A or G or T
                                                                       ret['H'] = ret['A'];
                                                                       ret['h'] = ret['A']; // A or C or T
                                                                       ret['V'] = ret['A'];
                                                                       ret['v'] = ret['A']; // A or C or G

                                                                       return ret;
                                                                   }()};
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::dna4.
 * \relates dna4
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna4_vector = std::vector<dna4>;

inline namespace literals
{

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Nucleotide literals
 * \{
 */

/*!\brief The seqan3::dna4 char literal.
 * \relatesalso seqan3::dna4
 * \returns seqan3::dna4
 * \details
 * \anchor seqan3_dna4_char_literal
 *
 * You can use this char literal to assign a seqan3::dna4 character:
 * \snippet test/snippet/alphabet/nucleotide/dna4_char_literal.cpp main
 *
 * \stableapi{Since version 3.1.}
 */
constexpr dna4 operator""_dna4(char const c) noexcept
{
    return dna4{}.assign_char(c);
}

/*!\brief The seqan3::dna4 string literal.
 * \relatesalso seqan3::dna4
 * \returns seqan3::dna4_vector
 * \anchor seqan3_dna4_string_literal
 *
 * You can use this string literal to easily assign to dna4_vector:
 * \include test/snippet/alphabet/nucleotide/dna4_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL dna4_vector operator""_dna4(char const * s, std::size_t n)
{
    dna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
