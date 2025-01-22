// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::rna4, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>

// ------------------------------------------------------------------
// rna4
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The four letter RNA alphabet of A,C,G,U.
 * \ingroup alphabet_nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * This alphabet has the same internal representation as seqan3::dna4, the only difference is that it prints 'U' on
 * character conversion instead of 'T'. You can assign between values of seqan3::dna4 and seqan3::rna4.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal \ref seqan3_rna4_char_literal "'A'_rna4" or use the
 * function seqan3::rna4::assign_char().
 *
 * \include test/snippet/alphabet/nucleotide/rna4.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class rna4 : public nucleotide_base<rna4, 4>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<rna4, 4>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond
    //!\brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr rna4() noexcept = default;                         //!< Defaulted.
    constexpr rna4(rna4 const &) noexcept = default;             //!< Defaulted.
    constexpr rna4(rna4 &&) noexcept = default;                  //!< Defaulted.
    constexpr rna4 & operator=(rna4 const &) noexcept = default; //!< Defaulted.
    constexpr rna4 & operator=(rna4 &&) noexcept = default;      //!< Defaulted.
    ~rna4() noexcept = default;                                  //!< Defaulted.

    using base_t::base_t;

    /*!\brief Allow implicit construction from seqan3::dna4 of the same size.
     * \details
     *
     * \include{doc} doc/fragments/rna4_implicit_conversion_from_dna4.md
     *
     * \stableapi{Since version 3.1.}
     */
    template <std::same_as<dna4> t>
    constexpr rna4(t const & r) noexcept
    {
        assign_rank(r.to_rank());
    }
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'C', 'G', 'U'};

    //!\copydoc seqan3::dna4::rank_complement
    static constexpr rank_type rank_complement(rank_type const rank)
    {
        return dna4::rank_complement(rank);
    }

    //!\copydoc seqan3::dna4::rank_to_char
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        return dna4::char_to_rank(chr);
    }
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::rna4.
 * \relates rna4
 * \details
 * \stableapi{Since version 3.1.}
 */
using rna4_vector = std::vector<rna4>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Nucleotide literals
 * \{
 */
/*!\brief The seqan3::rna4 char literal.
 * \relatesalso seqan3::rna4
 * \returns seqan3::rna4
 * \details
 * \anchor seqan3_rna4_char_literal
 *
 * You can use this char literal to assign a seqan3::rna4 character:
 * \snippet test/snippet/alphabet/nucleotide/rna4_char_literal.cpp main
 *
 * \stableapi{Since version 3.1.}
 */
constexpr rna4 operator""_rna4(char const c) noexcept
{
    return rna4{}.assign_char(c);
}

/*!\brief The seqan3::rna4 string literal.
 * \relatesalso seqan3::rna4
 * \returns seqan3::rna4_vector
 *
 * You can use this string literal to easily assign to rna4_vector:
 * \include test/snippet/alphabet/nucleotide/rna4_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL rna4_vector operator""_rna4(char const * s, std::size_t n)
{
    rna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
