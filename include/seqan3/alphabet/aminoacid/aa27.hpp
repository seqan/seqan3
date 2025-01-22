// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::aa27, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/aminoacid/aminoacid_base.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

namespace seqan3
{
/*!\brief The twenty-seven letter amino acid alphabet.
 * \ingroup alphabet_aminoacid
 * \implements seqan3::aminoacid_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * The alphabet consists of letters A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X,
 * Y, Z, *
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal \ref seqan3_aa27_char_literal "'X'_aa27" or use the
 * function seqan3::aa27::assign_char().
 *
 * \include test/snippet/alphabet/aminoacid/aa27.cpp
 *
 * \stableapi{Since version 3.1.}
 */

class aa27 : public aminoacid_base<aa27, 27>
{
private:
    //!\brief The base class.
    using base_t = aminoacid_base<aa27, 27>;

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
    constexpr aa27() noexcept = default;                         //!< Defaulted.
    constexpr aa27(aa27 const &) noexcept = default;             //!< Defaulted.
    constexpr aa27(aa27 &&) noexcept = default;                  //!< Defaulted.
    constexpr aa27 & operator=(aa27 const &) noexcept = default; //!< Defaulted.
    constexpr aa27 & operator=(aa27 &&) noexcept = default;      //!< Defaulted.
    ~aa27() noexcept = default;                                  //!< Defaulted.

    using base_t::base_t;
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
                                                                 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
                                                                 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '*'};

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
    static constexpr std::array<rank_type, 256> char_to_rank_table{
        []() constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with 'X' (UNKNOWN)
            ret.fill(23u);

            // reverse mapping for characters and their lowercase
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[static_cast<rank_type>(rank_to_char_table[rnk])] = rnk;
                ret[static_cast<rank_type>(to_lower(rank_to_char_table[rnk]))] = rnk;
            }

            return ret;
        }()};
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::aa27.
 * \relates aa27
 *
 * \stableapi{Since version 3.1.}
 */
using aa27_vector = std::vector<aa27>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Literals
 * \{
 */
/*!\brief The seqan3::aa27 char literal.
 * \param[in] c The character to assign.
 * \relates seqan3::aa27
 * \returns seqan3::aa27
 * \anchor seqan3_aa27_char_literal
 *
 * You can use this char literal to assign a seqan3::aa27 character:
 * \include test/snippet/alphabet/aminoacid/aa27_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr aa27 operator""_aa27(char const c) noexcept
{
    return aa27{}.assign_char(c);
}

/*!\brief The seqan3::aa27 string literal.
 * \param[in] s A pointer to the character string to assign.
 * \param[in] n The size of the character string to assign.
 * \relates seqan3::aa27
 * \returns seqan3::aa27_vector
 *
 * You can use this string literal to easily assign to aa27_vector:
 * \include test/snippet/alphabet/aminoacid/aa27_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL aa27_vector operator""_aa27(char const * const s, size_t const n)
{
    aa27_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
