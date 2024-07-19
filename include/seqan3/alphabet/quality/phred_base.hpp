// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::phred42 quality scores.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <algorithm>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/quality/concept.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that refines seqan3::alphabet_base and is used by the quality alphabets.
 * \ingroup alphabet_quality
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 * \details
 * \stableapi{Since version 3.1.}
 */
template <typename derived_type, size_t size>
class phred_base : public alphabet_base<derived_type, size, char>
{
public:
    /*!\name Member types
     * \{
     */
    /*!\brief The integer representation of the quality score.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using phred_type = int8_t;
    //!\}

private:
    //!\brief The base type.
    using base_t = alphabet_base<derived_type, size, char>;

    /*!\brief Befriend the base type so it can access seqan3::alphabet_base::char_to_rank
     *        and seqan3::alphabet_base::rank_to_char.
     */
    friend base_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr phred_base() noexcept = default;                               //!< Defaulted.
    constexpr phred_base(phred_base const &) noexcept = default;             //!< Defaulted.
    constexpr phred_base(phred_base &&) noexcept = default;                  //!< Defaulted.
    constexpr phred_base & operator=(phred_base const &) noexcept = default; //!< Defaulted.
    constexpr phred_base & operator=(phred_base &&) noexcept = default;      //!< Defaulted.
    ~phred_base() noexcept = default;                                        //!< Defaulted.

    //!\brief Befriend the derived_type so it can instantiate.
    friend derived_type;

public:
    // Import from base type:
    using base_t::alphabet_size;
    using base_t::assign_rank;
    using base_t::to_rank;
    using typename base_t::char_type;
    using typename base_t::rank_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    // This constructor needs to be public, because constructor templates are not inherited otherwise
    /*!\brief Allow explicit construction from any other quality type by means of the Phred score representation.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename other_qual_type>
        requires (!std::same_as<phred_base, other_qual_type>) && (!std::same_as<derived_type, other_qual_type>)
              && quality_alphabet<other_qual_type>
    explicit constexpr phred_base(other_qual_type const & other) noexcept
    {
        assign_phred_to(seqan3::to_phred(other), static_cast<derived_type &>(*this));
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the alphabet's value in Phred score representation.
     *
     * \see quality
     *
     * \stableapi{Since version 3.1.}
     */
    constexpr phred_type to_phred() const noexcept
    {
        return rank_to_phred[to_rank()];
    }
    //!\}

    /*!\name Write functions
     * \{
     */

    /*!\brief Assign from the numeric Phred score value.
     *
     * \details
     *
     * Satisfies the seqan3::writable_quality_alphabet requirement via the seqan3::assign_rank_to()
     * wrapper.
     *
     * \see quality
     *
     * ### Complexity
     *
     * Constant.
     *
     * \stableapi{Since version 3.1.}
     */
    constexpr derived_type & assign_phred(phred_type const p) noexcept
    {
        return assign_rank(phred_to_rank[static_cast<rank_type>(p)]);
    }
    //!\}

private:
    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        int64_t difference = static_cast<int64_t>(chr) - static_cast<int64_t>(derived_type::offset_char);
        return std::clamp<int64_t>(difference, 0, alphabet_size - 1);
    }

    //!\copydoc seqan3::dna4::rank_to_char
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank + derived_type::offset_char;
    }

    //!\brief Phred to rank conversion table.
    static constexpr std::array<rank_type, 256> phred_to_rank{
        []() constexpr
        {
            std::array<rank_type, 256> ret{};

            for (int64_t i = std::numeric_limits<phred_type>::lowest(); i <= std::numeric_limits<phred_type>::max();
                 ++i)
            {
                if (i < derived_type::offset_phred) // map too-small to smallest possible
                    ret[static_cast<rank_type>(i)] = 0;
                else if (i >= derived_type::offset_phred + alphabet_size) // map too-large to highest possible
                    ret[static_cast<rank_type>(i)] = alphabet_size - 1;
                else // map valid range to identity
                    ret[static_cast<rank_type>(i)] = i - derived_type::offset_phred;
            }

            return ret;
        }()};

    //!\brief Rank to phred conversion table.
    static constexpr std::array<phred_type, alphabet_size> rank_to_phred{
        []() constexpr
        {
            std::array<phred_type, alphabet_size> ret{};

            for (size_t i = 0; i < alphabet_size; ++i)
                ret[i] = i + derived_type::offset_phred;

            return ret;
        }()};
};

} // namespace seqan3
