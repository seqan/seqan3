// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Free function/metafunction wrappers for alphabets with member functions/types.
 *
 * This shall not need be included manually, just include `alphabet/concept.hpp`.
 */

#pragma once

#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/io/stream/char_operations.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that refines seqan3::alphabet_base and is used by the amino acids.
 * \ingroup aminoacid
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 */
template <typename derived_type, auto size>
class aminoacid_base : public alphabet_base<derived_type, size, char>
{
private:
    //!\brief Type of the base class.
    using base_t = alphabet_base<derived_type, size, char>;

    //!\brief Befriend the base class.
    friend base_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aminoacid_base()                                   noexcept = default; //!< Defaulted.
    constexpr aminoacid_base(aminoacid_base const &)             noexcept = default; //!< Defaulted.
    constexpr aminoacid_base(aminoacid_base &&)                  noexcept = default; //!< Defaulted.
    constexpr aminoacid_base & operator=(aminoacid_base const &) noexcept = default; //!< Defaulted.
    constexpr aminoacid_base & operator=(aminoacid_base &&)      noexcept = default; //!< Defaulted.
    ~aminoacid_base()                                            noexcept = default; //!< Defaulted.
    //!\}

    //!\brief Befriend the derived class so it can instantiate.
    friend derived_type;

public:

    // Import from base:
    using typename base_t::char_type;
    using typename base_t::rank_type;
    using base_t::value_size;
    using base_t::to_rank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    // This constructor needs to be public, because constructor templates are not inherited otherwise
    //!\brief Allow explicit construction from any other aminoacid type and convert via the character representation.
    template <typename other_aa_type>
    //!\cond
        requires !std::Same<aminoacid_base, other_aa_type> &&
                 !std::Same<derived_type, other_aa_type> &&
                 AminoacidAlphabet<other_aa_type>
    //!\endcond
    explicit constexpr aminoacid_base(other_aa_type const & other) noexcept
    {
        using seqan3::to_rank;
        static_cast<derived_type &>(*this) =
            detail::convert_through_char_representation<derived_type, other_aa_type>[to_rank(other)];
    }
    //!\}

    /*!\brief Validate whether a character value has a one-to-one mapping to an alphabet value.
     *
     * \details
     *
     * Models the seqan3::Semialphabet::char_is_valid_for() requirement via the seqan3::char_is_valid_for()
     * wrapper.
     *
     * Behaviour specific to amino acids: True also for lower case letters that silently convert to their upper case.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
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

            return ret;
        }()
    };
};

} // namespace seqan3
