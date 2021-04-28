// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::dna15, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dna15
// ------------------------------------------------------------------

namespace seqan3
{

class rna15;

/*!\brief The 15 letter DNA alphabet, containing all IUPAC smybols minus the gap.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * Note that you can assign 'U' as a character to dna15 and it will silently
 * be converted to 'T'.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal or use the
 * function seqan3::dna15::assign_char().
 *
 * \include test/snippet/alphabet/nucleotide/dna15.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class dna15 : public nucleotide_base<dna15, 15>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna15, 15>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond
    //!\brief Befriend seqan3::rna15 so it can copy #char_to_rank.
    friend rna15;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna15()                           noexcept = default; //!< Defaulted.
    constexpr dna15(dna15 const &)              noexcept = default; //!< Defaulted.
    constexpr dna15(dna15 &&)                   noexcept = default; //!< Defaulted.
    constexpr dna15 & operator=(dna15 const &)  noexcept = default; //!< Defaulted.
    constexpr dna15 & operator=(dna15 &&)       noexcept = default; //!< Defaulted.
    ~dna15()                                    noexcept = default; //!< Defaulted.

    using base_t::base_t;

    /*!\brief Allow implicit construction from dna/rna of the same size.
     * \details
     * \stableapi{Since version 3.1.}
     */
    template <std::same_as<rna15> t>    // Accept incomplete type
    constexpr dna15(t const & r) noexcept
    {
        assign_rank(r.to_rank());
    }
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]
    {
        'A',
        'B',
        'C',
        'D',
        'G',
        'H',
        'K',
        'M',
        'N',
        'R',
        'S',
        'T',
        'V',
        'W',
        'Y'
    };

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr std::array<rank_type, 256> char_to_rank_table
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = 8; // rank of 'N'

            // reverse mapping for characters and their lowercase
            for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[rank_to_char_table[rnk]] = rnk;
                ret[to_lower(rank_to_char_table[rnk])] = rnk;
            }

            // set U equal to T
            ret['U'] = ret['T']; ret['u'] = ret['t'];

            return ret;
        }()
    };

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

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<dna15, alphabet_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for an std::vector of seqan3::dna15.
 * \relates dna15
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna15_vector = std::vector<dna15>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Nucleotide literals
 * \{
 */
/*!\brief The seqan3::dna15 char literal.
 * \relatesalso seqan3::dna15
 * \returns seqan3::dna15
 * \details
 *
 * You can use this char literal to assign a seqan3::dna15 character:
 * \include test/snippet/alphabet/nucleotide/dna15_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr dna15 operator""_dna15(char const c) noexcept
{
    return dna15{}.assign_char(c);
}

/*!\brief The seqan3::dna15 string literal.
 * \relatesalso seqan3::dna15
 * \returns seqan3::dna15_vector
 *
 * You can use this string literal to easily assign to dna15_vector:
 * \include test/snippet/alphabet/nucleotide/dna15_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
inline dna15_vector operator""_dna15(char const * s, std::size_t n)
{
    dna15_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // inline namespace literals

// ------------------------------------------------------------------
// dna15 (deferred definition)
// ------------------------------------------------------------------

constexpr std::array<dna15, dna15::alphabet_size> dna15::complement_table
{
    'T'_dna15,    // complement of 'A'_dna15
    'V'_dna15,    // complement of 'B'_dna15
    'G'_dna15,    // complement of 'C'_dna15
    'H'_dna15,    // complement of 'D'_dna15
    'C'_dna15,    // complement of 'G'_dna15
    'D'_dna15,    // complement of 'H'_dna15
    'M'_dna15,    // complement of 'K'_dna15
    'K'_dna15,    // complement of 'M'_dna15
    'N'_dna15,    // complement of 'N'_dna15
    'Y'_dna15,    // complement of 'R'_dna15
    'S'_dna15,    // complement of 'S'_dna15
    'A'_dna15,    // complement of 'T'_dna15
    'B'_dna15,    // complement of 'V'_dna15
    'W'_dna15,    // complement of 'W'_dna15
    'R'_dna15     // complement of 'Y'_dna15
};

} // namespace seqan3
