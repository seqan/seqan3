// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Extends a given alphabet with the mask alphabet.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/mask/mask.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

namespace seqan3
{
/*!\brief Implementation of a masked composite, which extends a given alphabet with a mask.
 * \ingroup alphabet_mask
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 *
 * \tparam sequence_alphabet_t Type of the first letter; must satisfy seqan3::writable_alphabet and std::regular.
 *
 * \details
 *
 * The masked composite represents a seqan3::alphabet_tuple_base of any given alphabet with the masked alphabet. It
 * allows one to specify which portions of a sequence should be masked, without losing additional information by
 * replacing the sequence directly.
 *
 * \include test/snippet/alphabet/mask/masked.cpp
 *
 * \see \link mask Mask submodule \endlink, it contains an explanation of hard-masking (UNKNOWN character) and
 *      soft-masking (lower/upper case letters).
 *
 * \stableapi{Since version 3.1.}
 */
template <writable_alphabet sequence_alphabet_t>
    requires std::regular<sequence_alphabet_t>
class masked : public alphabet_tuple_base<masked<sequence_alphabet_t>, sequence_alphabet_t, mask>
{
private:
    //!\brief The base type.
    using base_t = alphabet_tuple_base<masked<sequence_alphabet_t>, sequence_alphabet_t, mask>;

public:
    /*!\brief First template parameter as member type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using sequence_alphabet_type = sequence_alphabet_t;

    /*!\brief Equals the char_type of sequence_alphabet_type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using char_type = alphabet_char_t<sequence_alphabet_type>;

    using base_t::alphabet_size;
    using typename base_t::rank_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr masked() = default;                           //!< Defaulted.
    constexpr masked(masked const &) = default;             //!< Defaulted.
    constexpr masked(masked &&) = default;                  //!< Defaulted.
    constexpr masked & operator=(masked const &) = default; //!< Defaulted.
    constexpr masked & operator=(masked &&) = default;      //!< Defaulted.
    ~masked() = default;                                    //!< Defaulted.

    using base_t::base_t; // Inherit non-default constructors
    //!\}

    // Inherit operators from base
    using base_t::operator=;

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from a character.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr masked & assign_char(char_type const c) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        base_t::assign_rank(char_to_rank[static_cast<index_t>(c)]);
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return a character.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr char_type to_char() const noexcept
    {
        return rank_to_char[base_t::to_rank()];
    }
    //!\}

protected:
    //!\brief Rank to char conversion table.
    static constexpr std::array<char_type, alphabet_size> rank_to_char{
        []() constexpr
        {
            std::array<char_type, alphabet_size> ret{};

            for (size_t i = 0; i < alphabet_size; ++i)
            {
                ret[i] = (i < alphabet_size / 2)
                           ? seqan3::to_char(seqan3::assign_rank_to(i, sequence_alphabet_type{}))
                           : to_lower(seqan3::to_char(seqan3::assign_rank_to(i / 2, sequence_alphabet_type{})));
            }

            return ret;
        }()};

    //!\brief Char to rank conversion table.
    static constexpr std::array<rank_type, detail::size_in_values_v<char_type>> char_to_rank{
        []() constexpr
        {
            std::array<rank_type, detail::size_in_values_v<char_type>> ret{};

            for (size_t i = 0; i < 256; ++i)
            {
                char_type const c = static_cast<char_type>(i);

                ret[i] = is_lower(c) ? seqan3::to_rank(seqan3::assign_char_to(c, sequence_alphabet_type{})) * 2
                                     : seqan3::to_rank(seqan3::assign_char_to(c, sequence_alphabet_type{}));
            }

            return ret;
        }()};
};

//!\brief Type deduction guide enables usage of masked without specifying template args.
//!\relates masked
template <typename sequence_alphabet_type>
masked(sequence_alphabet_type &&, mask const &) -> masked<std::decay_t<sequence_alphabet_type>>;

} //namespace seqan3
