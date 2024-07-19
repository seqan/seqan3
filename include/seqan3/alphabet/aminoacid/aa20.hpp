// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Provides seqan3::aa20, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/aminoacid/aminoacid_base.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

namespace seqan3
{

/*!\brief The canonical amino acid alphabet.
 * \ingroup alphabet_aminoacid
 * \implements seqan3::aminoacid_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * The alphabet consists of letters A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
 *
 * The alphabet may be brace initialized from the static letter members (see above). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 * *Note:* Letters which belong in the extended alphabet will be automatically converted based on the frequency
 * of their options.\n Terminator characters are converted to W, because the most commonly occurring
 * stop codon in higher eukaryotes is UGA<sup>2</sup>.
 * Anything unknown is converted to S, because it occurs most frequently across 53 vertebrates<sup>1</sup>.
 *
 * |Input Letter  |Converts to  |
 * |--------------|-------------|
 * |B             |D<sup>1</sup>|
 * |J             |L<sup>1</sup>|
 * |O             |L<sup>1</sup>|
 * |U             |C<sup>1</sup>|
 * |Z             |E<sup>1</sup>|
 * |X (Unknown)   |S<sup>1</sup>|
 * |* (Terminator)|W<sup>2</sup>|
 * <sup><b>1</b></sup>King, J. L., & Jukes, T. H. (1969). Non-Darwinian Evolution.
 * Science, 164(3881), 788-798. doi:10.1126/science.164.3881.788\n
 * <sup><b>2</b></sup>Trotta, E. (2016). Selective forces and mutational biases drive stop codon usage
 * in the human genome: a comparison with sense codon usage.
 * BMC Genomics, 17, 366. https://doi.org/10.1186/s12864-016-2692-4
 *
 * \include test/snippet/alphabet/aminoacid/aa20.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class aa20 : public aminoacid_base<aa20, 20>
{
private:
    //!\brief The base class.
    using base_t = aminoacid_base<aa20, 20>;

    //!\brief Befriend seqan3::aminoacid_base.
    friend base_t;
    //!\cond
    //!\brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aa20() noexcept = default;                         //!< Defaulted.
    constexpr aa20(aa20 const &) noexcept = default;             //!< Defaulted.
    constexpr aa20(aa20 &&) noexcept = default;                  //!< Defaulted.
    constexpr aa20 & operator=(aa20 const &) noexcept = default; //!< Defaulted.
    constexpr aa20 & operator=(aa20 &&) noexcept = default;      //!< Defaulted.
    ~aa20() noexcept = default;                                  //!< Defaulted.

    using base_t::base_t;
    //!\}

private:
    //!\copydoc seqan3::aa27::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                                                                 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};

    //!\copydoc seqan3::aa27::rank_to_char
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    //!\copydoc seqan3::aa27::char_to_rank
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
    }

    //!\copydoc seqan3::aa27::char_to_rank_table
    static constexpr std::array<rank_type, 256> char_to_rank_table{
        []() constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with 'S' because that appears most frequently
            ret.fill(15u);

            // reverse mapping for characters and their lowercase
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[static_cast<rank_type>(rank_to_char_table[rnk])] = rnk;
                ret[static_cast<rank_type>(to_lower(rank_to_char_table[rnk]))] = rnk;
            }

            ret['B'] = ret['D'];
            ret['b'] = ret['D']; // Convert b (either D/N) to D, since D occurs more frequently.
            ret['J'] = ret['L'];
            ret['j'] = ret['L']; // Convert j (either I/L) to L, since L occurs more frequently.
            ret['O'] = ret['L'];
            ret['o'] = ret['L']; // Convert Pyrrolysine to lysine.
            ret['U'] = ret['C'];
            ret['u'] = ret['C']; // Convert Selenocysteine to cysteine.
            ret['X'] = ret['S'];
            ret['x'] = ret['S']; // Convert unknown amino acids to serine.
            ret['Z'] = ret['E'];
            ret['z'] = ret['E']; // Convert z (either E/Q) to E, since E occurs more frequently.
            ret['*'] = ret['W']; // The most common stop codon is UGA. This is most similar to a Tryptophan.
            return ret;
        }()};
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::aa20.
 * \relates aa20
 *
 * \stableapi{Since version 3.1.}
 */
using aa20_vector = std::vector<aa20>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Literals
 * \{
 */
/*!\brief The seqan3::aa20 char literal.
 * \param[in] c The character to assign.
 * \relates seqan3::aa20
 * \returns seqan3::aa20
 *
 * You can use this char literal to assign a seqan3::aa20 character:
 * \include test/snippet/alphabet/aminoacid/aa20_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr aa20 operator""_aa20(char const c) noexcept
{
    return aa20{}.assign_char(c);
}

/*!\brief The seqan3::aa20  string literal.
 * \param[in] s A pointer to the character string to assign.
 * \param[in] n The size of the character string to assign.
 * \relates seqan3::aa20
 * \returns seqan3::aa20_vector
 *
 * You can use this string literal to easily assign to aa20_vector:
 * \include test/snippet/alphabet/aminoacid/aa20_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL aa20_vector operator""_aa20(char const * const s, size_t const n)
{
    aa20_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
