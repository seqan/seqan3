// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Introduces the seqan3::exposition_only::cigar_operation alphabet.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/alphabet_base.hpp>

// ------------------------------------------------------------------
// cigar_operation
// ------------------------------------------------------------------

namespace seqan3::exposition_only
{

/*!\brief The actual implementation of seqan3::cigar::operation for documentation purposes only.
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 *
 * \note This class only exists because of technical reasons,
 *       please always use seqan3::cigar::operation instead of seqan3::exposition_only::cigar_operation.
 *
 * \if DEV
 * \note We cannot declare seqan3::cigar::operation in-class, because we need to specify the second tuple element of
 *       seqan3::alphabet_tuple_base before we actually can declare it in-class. This is a trade-off to make
 *       seqan3::cigar a non-template class.
 * \endif
 *
 * \sa seqan3::exposition_only for an explanation on exposition-only.
 *
 * \noapi{Please always use seqan3::cigar::operation. The API-Stability of all members documented in this class applies
 *        directly to seqan3::cigar::operation.}
 */
class cigar_operation : public alphabet_base<cigar_operation, 9, char>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<cigar_operation, 9, char>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr cigar_operation() noexcept = default;                                    //!< Defaulted.
    constexpr cigar_operation(cigar_operation const &) noexcept = default;             //!< Defaulted.
    constexpr cigar_operation(cigar_operation &&) noexcept = default;                  //!< Defaulted.
    constexpr cigar_operation & operator=(cigar_operation const &) noexcept = default; //!< Defaulted.
    constexpr cigar_operation & operator=(cigar_operation &&) noexcept = default;      //!< Defaulted.
    ~cigar_operation() noexcept = default;                                             //!< Defaulted.

    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'M', 'D', 'I', 'S', 'H', 'N', 'P', 'X', '='};

    //!\copydoc seqan3::dna4::char_to_rank_table
    static constexpr std::array<rank_type, 256> char_to_rank_table{[]() constexpr
                                                                   {
                                                                       std::array<rank_type, 256> ret{};

                                                                       // reverse mapping for characters
                                                                       for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
                                                                       {
                                                                           ret[rank_to_char_table[rnk]] = rnk;
                                                                       }

                                                                       return ret;
                                                                   }()}; // namespace seqan3::exposition_only

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
};

} // namespace seqan3::exposition_only
