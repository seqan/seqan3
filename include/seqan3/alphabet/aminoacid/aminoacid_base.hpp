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

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface.
     */
    //!\{
    static derived_type constexpr A          = assign_char(derived_type{}, 'A');
    static derived_type constexpr B          = assign_char(derived_type{}, 'B');
    static derived_type constexpr C          = assign_char(derived_type{}, 'C');
    static derived_type constexpr D          = assign_char(derived_type{}, 'D');
    static derived_type constexpr E          = assign_char(derived_type{}, 'E');
    static derived_type constexpr F          = assign_char(derived_type{}, 'F');
    static derived_type constexpr G          = assign_char(derived_type{}, 'G');
    static derived_type constexpr H          = assign_char(derived_type{}, 'H');
    static derived_type constexpr I          = assign_char(derived_type{}, 'I');
    static derived_type constexpr J          = assign_char(derived_type{}, 'J');
    static derived_type constexpr K          = assign_char(derived_type{}, 'K');
    static derived_type constexpr L          = assign_char(derived_type{}, 'L');
    static derived_type constexpr M          = assign_char(derived_type{}, 'M');
    static derived_type constexpr N          = assign_char(derived_type{}, 'N');
    static derived_type constexpr O          = assign_char(derived_type{}, 'O');
    static derived_type constexpr P          = assign_char(derived_type{}, 'P');
    static derived_type constexpr Q          = assign_char(derived_type{}, 'Q');
    static derived_type constexpr R          = assign_char(derived_type{}, 'R');
    static derived_type constexpr S          = assign_char(derived_type{}, 'S');
    static derived_type constexpr T          = assign_char(derived_type{}, 'T');
    static derived_type constexpr U          = assign_char(derived_type{}, 'U');
    static derived_type constexpr V          = assign_char(derived_type{}, 'V');
    static derived_type constexpr W          = assign_char(derived_type{}, 'W');
    static derived_type constexpr X          = assign_char(derived_type{}, 'X');
    static derived_type constexpr Y          = assign_char(derived_type{}, 'Y');
    static derived_type constexpr Z          = assign_char(derived_type{}, 'Z');
    static derived_type constexpr TERMINATOR = assign_char(derived_type{}, '*');
    static derived_type constexpr UNKNOWN    = X;
    //!\}
};

} // namespace seqan3
