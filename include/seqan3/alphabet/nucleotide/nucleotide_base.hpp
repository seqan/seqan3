// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::nucleotide_base.
 */

#pragma once

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/core/char_operations/transform.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that refines seqan3::alphabet_base and is used by the nucleotides.
 * \ingroup nucleotide
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 *
 * \details
 *
 * You can use this class to define your own nucleotide alphabet, but types are not required to be based on it to model
 * seqan3::nucleotide_alphabet, it is purely a way to avoid code duplication.
 *
 * In addition to the requirements of seqan3::alphabet_base, the derived type needs to define the following static
 * member variable (can be private):
 *
 *   * `static std::array<THAT_TYPE, alphabet_size> complement_table` that defines for every possible rank value
 *     the corresponding complement.
 */
template <typename derived_type, auto size>
class nucleotide_base : public alphabet_base<derived_type, size, char>
{
private:
    //!\brief Type of the base class.
    using base_t = alphabet_base<derived_type, size, char>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr nucleotide_base()                                    noexcept = default; //!< Defaulted.
    constexpr nucleotide_base(nucleotide_base const &)             noexcept = default; //!< Defaulted.
    constexpr nucleotide_base(nucleotide_base &&)                  noexcept = default; //!< Defaulted.
    constexpr nucleotide_base & operator=(nucleotide_base const &) noexcept = default; //!< Defaulted.
    constexpr nucleotide_base & operator=(nucleotide_base &&)      noexcept = default; //!< Defaulted.
    ~nucleotide_base()                                             noexcept = default; //!< Defaulted.
    //!\}

    //! Befriend the derived_type so it can instantiate.
    friend derived_type;

protected:
    // Import from base:
    using typename base_t::char_type;
    using typename base_t::rank_type;

public:
    using base_t::alphabet_size;
    using base_t::to_rank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    // This constructor needs to be public, because constructor templates are not inherited otherwise
    //!\brief Allow explicit construction from any other nucleotide type and convert via the character representation.
    template <typename other_nucl_type>
    //!\cond
        requires !std::same_as<nucleotide_base, other_nucl_type> &&
                 !std::same_as<derived_type, other_nucl_type> &&
                 nucleotide_alphabet<other_nucl_type>
    //!\endcond
    explicit constexpr nucleotide_base(other_nucl_type const & other) noexcept
    {
        static_cast<derived_type &>(*this) =
            detail::convert_through_char_representation<derived_type, other_nucl_type>[seqan3::to_rank(other)];
    }
    //!\}

    /*!\name Read functions
     * \{
     */

    /*!\brief Return the complement of the letter.
     *
     * \details
     *
     * See \ref nucleotide for the actual values.
     *
     * Provides an implementation for seqan3::complement, required to model seqan3::nucleotide_alphabet.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr derived_type complement() const noexcept
    {
        return derived_type::complement_table[to_rank()];
    }
    //!\}

    /*!\brief Validate whether a character value has a one-to-one mapping to an alphabet value.
     *
     * \details
     *
     * Satisfies the seqan3::semialphabet::char_is_valid_for() requirement via the seqan3::char_is_valid_for()
     * wrapper.
     *
     * Behaviour specific to nucleotides: True also for lower case letters that silently convert to their upper case
     * **and** true also for U/T respectively, e.g. 'U' is a valid character for seqan3::dna4, because its informational
     * content is identical to 'T'.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Guaranteed not to throw.
     */
    static constexpr bool char_is_valid(char_type const c) noexcept
    {
        return valid_char_table[static_cast<uint8_t>(c)];
    }

private:
    //!\brief Implementation of #char_is_valid().
    static constexpr std::array<bool, 256> valid_char_table
    {
        [] () constexpr
        {
            // init with false
            std::array<bool, 256> ret{};

            // the original valid chars and their lower cases
            for (uint8_t c : derived_type::rank_to_char)
            {
                ret[         c ] = true;
                ret[to_lower(c)] = true;
            }

            // U and T shall be accepted for all
            ret['U'] = true;
            ret['T'] = true;
            ret['u'] = true;
            ret['t'] = true;

            return ret;
        }()
    };
};

} // namespace seqan3
