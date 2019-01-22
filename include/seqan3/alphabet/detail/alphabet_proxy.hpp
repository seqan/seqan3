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
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/metafunction/transformation_trait_or.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{
#if 0 // this is the alphabet_proxy I want, but GCC won't give me:
template <typename derived_type, typename alphabet_type>
class alphabet_proxy : public alphabet_type
{
public:
    using base_t = alphabet_type;

    using base_t::value_size;
    using typename base_t::rank_type;
    using char_type  = detail::transformation_trait_or_t<underlying_char<alphabet_type>, void>;
    using phred_type = detail::transformation_trait_or_t<underlying_phred<alphabet_type>, void>;

    using char_type_virtual  = std::conditional_t<std::Same<char_type, void>, char, char_type>;
    using phred_type_virtual = std::conditional_t<std::Same<phred_type, void>, int8_t, phred_type>;

    constexpr alphabet_proxy() : base_t{} {}
    constexpr alphabet_proxy(alphabet_proxy const &) = default;
    constexpr alphabet_proxy(alphabet_proxy &&) = default;
    constexpr alphabet_proxy & operator=(alphabet_proxy const &) = default;
    constexpr alphabet_proxy & operator=(alphabet_proxy &&) = default;
    ~alphabet_proxy() = default;

    constexpr alphabet_proxy(alphabet_type const a) :
        base_t{a}
    {}

    constexpr alphabet_proxy & operator=(alphabet_type const & c) noexcept
    {
        base_t::assign_rank(seqan3::to_rank(c));
        static_cast<derived_type &>(*this).on_update(); // <- this invokes the actual proxy behaviour!
        return *this;
    }

    template <typename indirect_assignable_type>
        requires std::Assignable<alphabet_type &, indirect_assignable_type>
    constexpr alphabet_proxy & operator=(indirect_assignable_type const & c) noexcept
    {
        alphabet_type a{};
        a = c;
        return operator=(a);
    }

    constexpr alphabet_proxy & assign_char(char_type_virtual const c) noexcept
        requires !std::Same<char_type, void>
    {
        alphabet_type tmp{};
        using seqan3::assign_char;
        assign_char(tmp, c);
        return operator=(tmp);
    }

    constexpr alphabet_proxy & assign_rank(underlying_rank_t<alphabet_type> const r) noexcept
    {
        alphabet_type tmp{};
        using seqan3::assign_rank;
        assign_rank(tmp, r);
        return operator=(tmp);
    }

    constexpr alphabet_proxy & assign_phred(phred_type_virtual const c) noexcept
        requires !std::Same<phred_type, void>
    {
        alphabet_type tmp{};
        using seqan3::assign_phred;
        assign_phred(tmp, c);
        return operator=(tmp);
    }
};
#endif

#if 1// this is the one that works for most things, but not all

/*!\brief A CRTP-base that eases the definition of proxy types returned in place of regular alphabets.
 * \tparam derived_type  The CRTP parameter type.
 * \tparam alphabet_type The type of the alphabet that this proxy emulates.
 *
 * \details
 *
 * Certain containers and other data structure hold alphabet values in a non-standard way so they can convert
 * to that alphabet when being accessed, but cannot return a reference to the held value. These data structures
 * may instead return a *proxy* to the held value which still allows changing it (and updating the underlying data
 * structure to reflect this).
 *
 * This CRTP base facilitates the definition of such proxies. Most users of SeqAn will not need to understand the
 * details.
 *
 * This class ensures that the proxy itself also models seqan3::semi_alphabet_concept, seqan3::alphabet_concept,
 * seqan3::quality_concept, seqan3::nucleotide_concept and/or seqan3::aminoacid_concept if the emulated type models
 * these. This makes sure that function templates which accept the original, also accept the proxy. An exception
 * are multi-layered compositions of alphabets where the proxy currently does not support access via `get`.
 *
 * ### Implementation notes
 *
 * The derived type needs to provide an `.on_update()` member function that performs the changes in the underlying
 * data structure.
 *
 * See seqan3::bitcompressed_vector or seqan3::cartesian_composition for examples of how this class is used.
 */
template <typename derived_type, typename alphabet_type>
class alphabet_proxy : public alphabet_base<derived_type,
                                            alphabet_size_v<alphabet_type>,
                                            detail::transformation_trait_or_t<underlying_char<alphabet_type>, void>>
{
private:
    //!\brief Type of the base class.
    using base_t =  alphabet_base<derived_type,
                                  alphabet_size_v<alphabet_type>,
                                  detail::transformation_trait_or_t<underlying_char<alphabet_type>, void>>;

    //!\brief Befriend the base type.
    friend base_t;

public:
    // Import from base:
    using base_t::value_size;
    using base_t::to_rank;

    /*!\name Member types
     * \{
     */
    using rank_type  = underlying_rank_t<alphabet_type>;
    using char_type  = detail::transformation_trait_or_t<underlying_char<alphabet_type>, void>;
    using phred_type = detail::transformation_trait_or_t<underlying_phred<alphabet_type>, void>;
    //!\}

private:
    //!\brief Never used, but required for valid definitions.
    using char_type_virtual  = std::conditional_t<std::Same<char_type, void>, char, char_type>;
    //!\brief Never used, but required for valid definitions.
    using phred_type_virtual = std::conditional_t<std::Same<phred_type, void>, int8_t, phred_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_proxy() noexcept : base_t{} {}
    constexpr alphabet_proxy(alphabet_proxy const &) = default;
    constexpr alphabet_proxy(alphabet_proxy &&) = default;
    constexpr alphabet_proxy & operator=(alphabet_proxy const &) = default;
    constexpr alphabet_proxy & operator=(alphabet_proxy &&) = default;
    ~alphabet_proxy() = default;

    //!\brief Construction from the emulated type.
    constexpr alphabet_proxy(alphabet_type const a) noexcept
    {
        base_t::assign_rank(seqan3::to_rank(a));
    }

    //!\brief Assigment from the emulated type. This function triggers the specialisation in the derived_type.
    constexpr derived_type & operator=(alphabet_type const & c) noexcept
    {
        base_t::assign_rank(seqan3::to_rank(c));
        static_cast<derived_type &>(*this).on_update(); // <- this invokes the actual proxy behaviour!
        return static_cast<derived_type &>(*this);
    }

    //!\brief Assignment from any type that the emulated type is assignable from.
    template <typename indirect_assignable_type>
    constexpr derived_type & operator=(indirect_assignable_type const & c) noexcept
        requires weakly_assignable_concept<alphabet_type, indirect_assignable_type>
    {
        alphabet_type a{};
        a = c;
        return operator=(a);
    }
    //!\}

    //!\brief Befriend the derived type so it can instantiate.
    friend derived_type;

public:
    /*!\name Write functions
     * \brief All of these call the emulated type's write functions and then delegate to
     *        the assignment operator which invokes derived behaviour.
     * \{
     */
    constexpr derived_type & assign_rank(underlying_rank_t<alphabet_type> const r) noexcept
    {
        alphabet_type tmp{};
        using seqan3::assign_rank;
        assign_rank(tmp, r);
        return operator=(tmp);
    }

    constexpr derived_type & assign_char(char_type_virtual const c) noexcept
        requires alphabet_concept<alphabet_type>
    {
        alphabet_type tmp{};
        using seqan3::assign_char;
        assign_char(tmp, c);
        return operator=(tmp);
    }

    derived_type & assign_char_strict(char_type_virtual const c)
        requires alphabet_concept<alphabet_type>
    {
        alphabet_type tmp{};
        using seqan3::assign_char_strict;
        assign_char_strict(tmp, c);
        return operator=(tmp);
    }

    constexpr derived_type & assign_phred(phred_type_virtual const c) noexcept
        requires quality_concept<alphabet_type>
    {
        alphabet_type tmp{};
        using seqan3::assign_phred;
        assign_phred(tmp, c);
        return operator=(tmp);
    }
    //!\}

    /*!\name Read functions
     * \brief All of these call the emulated type's read functions.
     * \{
     */
    //!\brief Implicit conversion to the emulated type.
    constexpr operator alphabet_type() const noexcept
    {
        using seqan3::assign_rank;
        return assign_rank(alphabet_type{}, to_rank());
    }

    constexpr char_type to_char() const noexcept
        requires alphabet_concept<alphabet_type>
    {
        using seqan3::to_char;
        return to_char(static_cast<alphabet_type>(*this));
    }

    constexpr phred_type to_phred() const noexcept
        requires quality_concept<alphabet_type>
    {
        using seqan3::to_phred;
        return to_phred(static_cast<alphabet_type>(*this));
    }

#if 0 // this currently causes GCC ICE in cartesian_composition test
    constexpr alphabet_type complement() const noexcept
        requires nucleotide_concept<alphabet_type>
    {
        using seqan3::complement;
        return complement(static_cast<alphabet_type>(*this));
    }
#endif

    //!\brief Delegate to the emulated type's validator.
    static constexpr bool char_is_valid(char_type_virtual const c) noexcept
    {
        using seqan3::char_is_valid_for;
        return char_is_valid_for<alphabet_type>(c);
    }
    //!\}
};

#endif
} // namespace seqan3
