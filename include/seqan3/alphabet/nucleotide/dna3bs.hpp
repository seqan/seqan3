// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::dna3bs, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dna3bs
// ------------------------------------------------------------------

namespace seqan3
{
/*!\brief The three letter reduced DNA alphabet for bisulfite sequencing mode (A,G,T(=C)).
 * \ingroup alphabet_nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * This alphabet represents a reduced version that can be used when dealing with bisulfite-converted data.
 * All 'C's are converted to a 'T' in order to allow comparison of normal sequences with
 * bisulfite-converted sequences.
 * For completeness, this nucleotide alphabet has a complement table, however, it is not recommended to use
 * it when dealing with bisulfite data because the complement of T is ambiguous in reads from bisulfite
 * sequencing. A 'T' can represent a true thymidine or an unmethylated 'C' that was converted into a 'T'.
 * Therefore, complementing a seqan3::dna3bs sequence will further reduce the alphabet to only 'T' and 'A', thereby
 * losing all information about 'G'. When working with bisulfite data, we recommend to create the reverse
 * complement of the seqan3::dna4 / seqan3::dna5 / seqan3::dna15 range first and convert to seqan3::dna3bs later. This
 * avoids simplifying the data by automatically setting 'A' as the complement of 'C'. As an example: The sequence
 * 'ACGTGC' in seqan3::dna4 would be 'ATGTGT' in seqan3::dna3bs. The complement of this seqan3::dna3bs sequence would
 * be 'TATATA', however when complementing the seqan3::dna4 sequence first and afterwards transforming it into
 * seqan3::dna3bs, it would be 'TGTATG' which preserves more information from the original sequence.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal \ref seqan3_dna3bs_char_literal "'A'_dna3bs" or use the
 * function seqan3::dna3bs::assign_char().
 *
 * \include test/snippet/alphabet/nucleotide/dna3bs.cpp
 *
 * \sa https://en.wikipedia.org/wiki/Bisulfite_sequencing
 *
 * \stableapi{Since version 3.1.}
 */
class dna3bs : public nucleotide_base<dna3bs, 3>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna3bs, 3>;

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
    constexpr dna3bs() noexcept = default;                           //!< Defaulted.
    constexpr dna3bs(dna3bs const &) noexcept = default;             //!< Defaulted.
    constexpr dna3bs(dna3bs &&) noexcept = default;                  //!< Defaulted.
    constexpr dna3bs & operator=(dna3bs const &) noexcept = default; //!< Defaulted.
    constexpr dna3bs & operator=(dna3bs &&) noexcept = default;      //!< Defaulted.
    ~dna3bs() noexcept = default;                                    //!< Defaulted.

    using base_t::base_t;
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'G', 'T'};

    //!\copydoc seqan3::dna4::rank_complement_table
    static constexpr rank_type rank_complement_table[alphabet_size]{
        2, // T is complement of 'A'_dna3bs
        2, // T is complement of 'G'_dna3bs
        0  // A is complement of 'T'_dna3bs
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

                                                                       // reverse mapping for characters and their lowercase
                                                                       for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
                                                                       {
                                                                           ret[rank_to_char_table[rnk]] = rnk;
                                                                           ret[to_lower(rank_to_char_table[rnk])] = rnk;
                                                                       }

                                                                       // set C and U equal to T
                                                                       ret['C'] = ret['T'];
                                                                       ret['c'] = ret['t'];
                                                                       ret['U'] = ret['T'];
                                                                       ret['u'] = ret['t'];

                                                                       // iupac characters get special treatment, because there is no N
                                                                       ret['R'] = ret['A'];
                                                                       ret['r'] = ret['A']; // A or G becomes A
                                                                       ret['Y'] = ret['T'];
                                                                       ret['y'] = ret['T']; // C or T becomes T
                                                                       ret['S'] = ret['T'];
                                                                       ret['s'] = ret['T']; // C or G becomes T
                                                                       ret['W'] = ret['A'];
                                                                       ret['w'] = ret['A']; // A or T becomes A
                                                                       ret['K'] = ret['G'];
                                                                       ret['k'] = ret['G']; // G or T becomes G
                                                                       ret['M'] = ret['A'];
                                                                       ret['m'] = ret['A']; // A or C becomes A
                                                                       ret['B'] = ret['T'];
                                                                       ret['b'] = ret['T']; // C or G or T becomes T
                                                                       ret['D'] = ret['A'];
                                                                       ret['d'] = ret['A']; // A or G or T becomes A
                                                                       ret['H'] = ret['A'];
                                                                       ret['h'] = ret['A']; // A or C or T becomes A
                                                                       ret['V'] = ret['A'];
                                                                       ret['v'] = ret['A']; // A or C or G  becomes A

                                                                       return ret;
                                                                   }()};
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::dna3bs.
 * \relates dna3bs
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna3bs_vector = std::vector<dna3bs>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Nucleotide literals
 * \{
 */
/*!\brief The seqan3::dna3bs char literal.
 * \relatesalso seqan3::dna3bs
 * \returns seqan3::dna3bs
 * \details
 * \anchor seqan3_dna3bs_char_literal
 *
 * You can use this char literal to assign a seqan3::dna3bs character:
 * \include test/snippet/alphabet/nucleotide/dna3bs_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr dna3bs operator""_dna3bs(char const c) noexcept
{
    return dna3bs{}.assign_char(c);
}

/*!\brief The seqan3::dna3bs string literal.
 * \relatesalso seqan3::dna3bs
 * \returns seqan3::dna3bs_vector
 *
 * You can use this string literal to easily assign to dna3bs_vector:
 * \include test/snippet/alphabet/nucleotide/dna3bs_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL dna3bs_vector operator""_dna3bs(char const * s, std::size_t n)
{
    dna3bs_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
