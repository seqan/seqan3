// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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
    constexpr aminoacid_base() : base_t{} {}
    constexpr aminoacid_base(aminoacid_base const &) = default;
    constexpr aminoacid_base(aminoacid_base &&) = default;
    constexpr aminoacid_base & operator=(aminoacid_base const &) = default;
    constexpr aminoacid_base & operator=(aminoacid_base &&) = default;
    ~aminoacid_base() = default;
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
                 aminoacid_concept<other_aa_type>
    //!\endcond
    explicit constexpr aminoacid_base(other_aa_type const & other) noexcept
    {
        using seqan3::to_rank;
        static_cast<derived_type &>(*this) =
            detail::convert_through_char_representation<derived_type, other_aa_type>[to_rank(other)];
    }
    //!\}
};

} // namespace seqan3
