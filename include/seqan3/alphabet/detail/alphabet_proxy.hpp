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
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/metafunction/transformation_trait_or.hpp>
#include <seqan3/std/concepts>

#if 0 // this is the alphabet_proxy I want, but GCC won't give me:
template <typename derived_type, typename alphabet_type>
class alphabet_proxy : public alphabet_type
{
public:
    using base_t = alphabet_type;

    using base_t::alphabet_size;
    using typename base_t::rank_type;
    using char_type  = detail::transformation_trait_or_t<std::type_identity<alphabet_char_t<<alphabet_type>, void>;
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
        using seqan3::assign_char_to;
        assign_char_to(c, tmp);
        return operator=(tmp);
    }

    constexpr alphabet_proxy & assign_rank(alphabet_rank_t<alphabet_type> const r) noexcept
    {
        alphabet_type tmp{};
        using seqan3::assign_rank_to;
        assign_rank_to(r, tmp);
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

namespace seqan3::detail
{

//!\brief Auxiliary transformation trait that is void for types by default.
//!\ingroup alphabet
template <typename t>
struct alphabet_char_type_or_void
{
    //!\brief The "return" type.
    using type = void;
};

//!\brief Auxiliary transformation trait specialisation that returns seqan3::alphabet_char_t.
//!\ingroup alphabet
template <Alphabet t>
struct alphabet_char_type_or_void<t>
{
    //!\brief The "return" type.
    using type = alphabet_char_t<t>;
};

//!\brief Auxiliary transformation trait shortcut.
//!\relates seqan3::alphabet_char_type_or_void
template <typename t>
using alphabet_char_type_or_void_t = typename alphabet_char_type_or_void<t>::type;

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief A CRTP-base that eases the definition of proxy types returned in place of regular alphabets.
 * \tparam derived_type  The CRTP parameter type.
 * \tparam alphabet_type The type of the alphabet that this proxy emulates.
 * \ingroup alphabet
 *
 * \details
 *
 * \noapi
 *
 * Certain containers and other data structure hold alphabet values in a non-standard way so they can convert
 * to that alphabet when being accessed, but cannot return a reference to the held value. These data structures
 * may instead return a *proxy* to the held value which still allows changing it (and updating the underlying data
 * structure to reflect this).
 *
 * This CRTP base facilitates the definition of such proxies. Most users of SeqAn will not need to understand the
 * details.
 *
 * This class ensures that the proxy itself also models seqan3::Semialphabet, seqan3::Alphabet,
 * seqan3::QualityAlphabet, seqan3::NucleotideAlphabet and/or seqan3::AminoacidAlphabet if the emulated type models
 * these. This makes sure that function templates which accept the original, also accept the proxy. An exception
 * are multi-layered composites of alphabets where the proxy currently does not support access via `get`.
 *
 * ### Implementation notes
 *
 * The derived type needs to provide an `.on_update()` member function that performs the changes in the underlying
 * data structure.
 *
 * See seqan3::bitcompressed_vector or seqan3::alphabet_tuple_base for examples of how this class is used.
 */
template <typename derived_type, WritableSemialphabet alphabet_type>
class alphabet_proxy : public alphabet_base<derived_type,
                                            alphabet_size<alphabet_type>,
                                            detail::alphabet_char_type_or_void_t<alphabet_type>>
{
private:
    //!\brief Type of the base class.
    using base_t =  alphabet_base<derived_type,
                                  alphabet_size<alphabet_type>,
                                  detail::alphabet_char_type_or_void_t<alphabet_type>>;

    //!\brief Befriend the base type.
    friend base_t;

public:
    // Import from base:
    using base_t::alphabet_size;
    using base_t::to_rank;

    /*!\name Member types
     * \{
     */
    using rank_type  = alphabet_rank_t<alphabet_type>;
    using char_type  = detail::alphabet_char_type_or_void_t<alphabet_type>;
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
    constexpr alphabet_proxy() noexcept : base_t{} {}                       //!< Defaulted.
    constexpr alphabet_proxy(alphabet_proxy const &) = default;             //!< Defaulted.
    constexpr alphabet_proxy(alphabet_proxy &&) = default;                  //!< Defaulted.
    constexpr alphabet_proxy & operator=(alphabet_proxy const &) = default; //!< Defaulted.
    constexpr alphabet_proxy & operator=(alphabet_proxy &&) = default;      //!< Defaulted.
    ~alphabet_proxy() = default;                                            //!< Defaulted.

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
        requires WeaklyAssignable<alphabet_type, indirect_assignable_type>
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
    constexpr derived_type & assign_rank(alphabet_rank_t<alphabet_type> const r) noexcept
    {
        alphabet_type tmp{};
        assign_rank_to(r, tmp);
        return operator=(tmp);
    }

    constexpr derived_type & assign_char(char_type_virtual const c) noexcept
        requires WritableAlphabet<alphabet_type>
    {
        alphabet_type tmp{};
        assign_char_to(c, tmp);
        return operator=(tmp);
    }

    constexpr derived_type & assign_phred(phred_type_virtual const c) noexcept
        requires QualityAlphabet<alphabet_type>
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
        return assign_rank_to(to_rank(), alphabet_type{});
    }

    constexpr char_type to_char() const noexcept
        requires Alphabet<alphabet_type>
    {
        /* (smehringer) Explicit conversion instead of static_cast:
         * See explanation in to_phred().
         */
        return seqan3::to_char(operator alphabet_type());
    }

    constexpr phred_type to_phred() const noexcept
        requires QualityAlphabet<alphabet_type>
    {
        using seqan3::to_phred;
        /* (smehringer) Explicit conversion instead of static_cast:
         * The tuple composite qualified returns a component_proxy which inherits from alphabet_proxy_base.
         * The qualified alphabet itself inherits from quality_base.
         * Now when accessing get<1>(seq_qual_alph) we want to call to_phred at some point because we want the quality,
         * therefore the to_phred function from alphabet_proxy is called, but this function did a static_cast to the
         * derived type which is calling the constructor from quality_base. Unfortunately now, the generic quality_base
         *  constructor uses assign_phred(static_cast<derived_type &>(*this), to_phred(other)); (here) which again
         * tries to call to_phred of the alphabet_proxy => infinite loop :boom:
         */
        return to_phred(operator alphabet_type());
    }

#if 0 // this currently causes GCC ICE in alphabet_tuple_base test
    constexpr alphabet_type complement() const noexcept
        requires NucleotideAlphabet<alphabet_type>
    {
        using seqan3::complement;
        return complement(static_cast<alphabet_type>(*this));
    }
#endif

    //!\brief Delegate to the emulated type's validator.
    static constexpr bool char_is_valid(char_type_virtual const c) noexcept
        requires WritableAlphabet<alphabet_type>
    {
        return char_is_valid_for<alphabet_type>(c);
    }
    //!\}
};

#endif
} // namespace seqan3
