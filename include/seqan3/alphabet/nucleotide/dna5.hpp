// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::dna5, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dna5
// ------------------------------------------------------------------

namespace seqan3
{

class rna5;

/*!\brief The five letter DNA alphabet of A,C,G,T and the unknown character N.
 * \ingroup alphabet_nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * Note that you can assign 'U' as a character to dna5 and it will silently
 * be converted to 'T'.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal \ref seqan3_dna5_char_literal "'A'_dna5" or use the
 * function seqan3::dna5::assign_char().
 *
 * \include test/snippet/alphabet/nucleotide/dna5.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class dna5 : public nucleotide_base<dna5, 5>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna5, 5>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond
    //!\brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond
    //!\brief Befriend seqan3::rna5 so it can copy #char_to_rank.
    friend rna5;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna5() noexcept = default;                         //!< Defaulted.
    constexpr dna5(dna5 const &) noexcept = default;             //!< Defaulted.
    constexpr dna5(dna5 &&) noexcept = default;                  //!< Defaulted.
    constexpr dna5 & operator=(dna5 const &) noexcept = default; //!< Defaulted.
    constexpr dna5 & operator=(dna5 &&) noexcept = default;      //!< Defaulted.
    ~dna5() noexcept = default;                                  //!< Defaulted.

    using base_t::base_t;

    /*!\brief Allow implicit construction from seqan3::rna5 of the same size.
     * \details
     *
     * \include{doc} doc/fragments/dna5_implicit_conversion_from_rna5.md
     *
     * \stableapi{Since version 3.1.}
     */
    template <std::same_as<rna5> t> // Accept incomplete type
    constexpr dna5(t const & r) noexcept
    {
        assign_rank(r.to_rank());
    }
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'C', 'G', 'N', 'T'};

    //!\copydoc seqan3::dna4::rank_complement_table
    static constexpr rank_type rank_complement_table[alphabet_size]{
        4, // T is complement of 'A'_dna5
        2, // G is complement of 'C'_dna5
        1, // C is complement of 'G'_dna5
        3, // N is complement of 'N'_dna5
        0  // A is complement of 'T'_dna5
    };

    //!\copydoc seqan3::dna4::rank_complement
    static constexpr rank_type rank_complement(rank_type const rank)
    {
        return rank_complement_table[rank];
    }

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
                                                                       std::array<rank_type, 256> ret{};

                                                                       ret.fill(3u); // initialize with UNKNOWN ('N')

                                                                       // reverse mapping for characters and their lowercase
                                                                       for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
                                                                       {
                                                                           ret[rank_to_char_table[rnk]] = rnk;
                                                                           ret[to_lower(rank_to_char_table[rnk])] = rnk;
                                                                       }

                                                                       // set U equal to T
                                                                       ret['U'] = ret['T'];
                                                                       ret['u'] = ret['t'];

                                                                       // iupac characters are implicitly "UNKNOWN"
                                                                       return ret;
                                                                   }()};
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::dna5.
 * \relates dna5
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna5_vector = std::vector<dna5>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Nucleotide literals
 * \{
 */
/*!\brief The seqan3::dna5 char literal.
 * \relatesalso seqan3::dna5
 * \returns seqan3::dna5
 * \details
 * \anchor seqan3_dna5_char_literal
 *
 * You can use this char literal to assign a seqan3::dna4 character:
 * \snippet test/snippet/alphabet/nucleotide/dna4_char_literal.cpp main
 *
 * \stableapi{Since version 3.1.}
 */
constexpr dna5 operator""_dna5(char const c) noexcept
{
    return dna5{}.assign_char(c);
}

/*!\brief The seqan3::dna5 string literal.
 * \relatesalso seqan3::dna5
 * \returns seqan3::dna5_vector
 *
 * You can use this string literal to easily assign to dna5_vector:
 * \include test/snippet/alphabet/nucleotide/dna5_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL dna5_vector operator""_dna5(char const * s, std::size_t n)
{
    dna5_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
