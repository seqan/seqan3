// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::rna5, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>

// ------------------------------------------------------------------
// rna5
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The five letter RNA alphabet of A,C,G,U and the unknown character N.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * This alphabet has the same internal representation as seqan3::dna5, the only difference is that it prints 'U' on
 * character conversion instead of 'T'. You can assign between values of seqan3::dna5 and seqan3::rna5.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal or use the
 * function seqan3::rna5::assign_char().
 *
 * \include test/snippet/alphabet/nucleotide/rna5.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class rna5 : public nucleotide_base<rna5, 5>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<rna5, 5>;

    //!\brief Befriend nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr rna5()                          noexcept = default; //!< Defaulted.
    constexpr rna5(rna5 const &)              noexcept = default; //!< Defaulted.
    constexpr rna5(rna5 &&)                   noexcept = default; //!< Defaulted.
    constexpr rna5 & operator=(rna5 const &)  noexcept = default; //!< Defaulted.
    constexpr rna5 & operator=(rna5 &&)       noexcept = default; //!< Defaulted.
    ~rna5()                                   noexcept = default; //!< Defaulted.

    using base_t::base_t;

    /*!\brief Allow implicit construction from dna/rna of the same size.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr rna5(dna5 const & r) noexcept
#if SEQAN3_WORKAROUND_GCC_90897
        requires true
#endif
    {
        assign_rank(r.to_rank());
    }
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]
    {
        'A',
        'C',
        'G',
        'N',
        'U'
    };

    //!\copydoc seqan3::dna4::rank_complement
    static constexpr rank_type rank_complement(rank_type const rank)
    {
        return dna5::rank_complement(rank);
    }

    //!\copydoc seqan3::dna4::rank_to_char
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        return dna5::char_to_rank(chr);
    }
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for an std::vector of seqan3::rna5.
 * \relates rna5
 * \details
 * \stableapi{Since version 3.1.}
 */
using rna5_vector = std::vector<rna5>;


// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Nucleotide literals
 * \{
 */
/*!\brief The seqan3::rna5 char literal.
 * \relatesalso seqan3::rna5
 * \returns seqan3::rna5
 * \details
 *
 * You can use this char literal to assign a seqan3::rna5 character:
 * \include test/snippet/alphabet/nucleotide/rna5_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr rna5 operator""_rna5(char const c) noexcept
{
    return rna5{}.assign_char(c);
}

/*!\brief The seqan3::rna5 string literal.
 * \relatesalso seqan3::rna5
 * \returns seqan3::rna5_vector
 *
 * You can use this string literal to easily assign to rna5_vector:
 * \include test/snippet/alphabet/nucleotide/rna5_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
inline rna5_vector operator""_rna5(char const * s, std::size_t n)
{
    rna5_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // inline namespace literals

} // namespace seqan3
