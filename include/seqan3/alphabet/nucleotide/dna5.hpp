// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::dna5, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/io/stream/char_operations.hpp>

// ------------------------------------------------------------------
// dna5
// ------------------------------------------------------------------

namespace seqan3
{

class rna5;

/*!\brief The five letter DNA alphabet of A,C,G,T and the unknown character N.
 * \ingroup nucleotide
 * \implements seqan3::NucleotideAlphabet
 * \implements seqan3::WritableAlphabet
 * \if DEV \implements seqan3::detail::WritableConstexprAlphabet \endif
 * \implements seqan3::TriviallyCopyable
 * \implements seqan3::StandardLayout
 * \implements std::Regular
 *
 * \details
 * Note that you can assign 'U' as a character to dna5 and it will silently
 * be converted to 'T'.
 *
 * The alphabet may be brace initialized from the static letter members. Note that you cannot
 * assign the alphabet by using letters of type `char`, but you instead have to use the
 * function seqan3::dna5::assign_char().
 *
 *\snippet test/snippet/alphabet/nucleotide/dna5.cpp code
 */
class dna5 : public nucleotide_base<dna5, 5>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna5, 5>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond
    //!\brief Befriend seqan3::rna5 so it can copy #char_to_rank.
    friend rna5;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna5()                          noexcept = default; //!< Defaulted.
    constexpr dna5(dna5 const &)              noexcept = default; //!< Defaulted.
    constexpr dna5(dna5 &&)                   noexcept = default; //!< Defaulted.
    constexpr dna5 & operator=(dna5 const &)  noexcept = default; //!< Defaulted.
    constexpr dna5 & operator=(dna5 &&)       noexcept = default; //!< Defaulted.
    ~dna5()                                   noexcept = default; //!< Defaulted.

    using base_t::base_t;

    //!\brief Allow implicit construction from dna/rna of the same size.
    template <std::Same<rna5> t>    // Accept incomplete type
    constexpr dna5(t const & r) noexcept
    {
        assign_rank(r.to_rank());
    }
    //!\}

protected:
    //!\privatesection

    //!\copydoc seqan3::dna4::rank_to_char
    static constexpr char_type rank_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'N',
        'T'
    };

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = 3; // == 'N'

            // reverse mapping for characters and their lowercase
            for (size_t rnk = 0u; rnk < value_size; ++rnk)
            {
                ret[         rank_to_char[rnk] ] = rnk;
                ret[to_lower(rank_to_char[rnk])] = rnk;
            }

            // set U equal to T
            ret['U'] = ret['T']; ret['u'] = ret['t'];

            // iupac characters are implicitly "UNKNOWN"
            return ret;
        }()
    };

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<dna5, value_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::dna5.
//!\relates dna5
using dna5_vector = std::vector<dna5>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::dna5 char literal.
 * \relates seqan3::dna5
 * \returns seqan3::dna5
 */
constexpr dna5 operator""_dna5(char const c) noexcept
{
    return dna5{}.assign_char(c);
}

/*!\brief The seqan3::dna5 string literal.
 * \relates seqan3::dna5
 * \returns seqan3::dna5_vector
 *
 * You can use this string literal to easily assign to dna5_vector:
 *
 * \snippet test/snippet/alphabet/nucleotide/dna5.cpp operator""_dna5
 *
 */
inline dna5_vector operator""_dna5(char const * s, std::size_t n)
{
    dna5_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// dna5 (deferred definition)
// ------------------------------------------------------------------

constexpr std::array<dna5, dna5::value_size> dna5::complement_table
{
    'T'_dna5,    // complement of 'A'_dna5
    'G'_dna5,    // complement of 'C'_dna5
    'C'_dna5,    // complement of 'G'_dna5
    'N'_dna5,    // complement of 'N'_dna5
    'A'_dna5     // complement of 'T'_dna5
};

} // namespace seqan3
