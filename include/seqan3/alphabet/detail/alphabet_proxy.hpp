// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_proxy.
 */

#pragma once

#include <concepts>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/detail/concept.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that eases the definition of proxy types returned in place of regular alphabets.
 * \tparam derived_type  The CRTP parameter type.
 * \tparam alphabet_type The type of the alphabet that this proxy emulates; must model at least
 *                       seqan3::writable_semialphabet and std::regular.
 * \ingroup alphabet
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
 * This class ensures that the proxy itself also models seqan3::semialphabet, seqan3::alphabet,
 * seqan3::quality_alphabet, seqan3::nucleotide_alphabet and/or seqan3::aminoacid_alphabet if the emulated type models
 * these. This makes sure that function templates which accept the original, also accept the proxy.
 *
 * ### Implementation notes
 *
 * The derived type needs to provide an `.on_update()` member function that performs the changes in the underlying
 * data structure.
 *
 * See seqan3::bitpacked_sequence or seqan3::alphabet_tuple_base for examples of how this class is used.
 *
 * \noapi{Exposition only}
 */
template <typename derived_type, writable_semialphabet alphabet_type>
    requires std::regular<alphabet_type>
class alphabet_proxy :
    public std::conditional_t<std::is_class_v<alphabet_type>,
                              alphabet_type,
                              alphabet_base<derived_type,
                                            alphabet_size<alphabet_type>,
                                            detail::valid_template_spec_or_t<void, alphabet_char_t, alphabet_type>>>
{
private:
    //!\brief Type of the base class.
    using base_t =
        std::conditional_t<std::is_class_v<alphabet_type>,
                           alphabet_type,              // inherit from emulated type if possible
                           alphabet_base<derived_type, // else: alphabet_base
                                         alphabet_size<alphabet_type>,
                                         detail::valid_template_spec_or_t<void, alphabet_char_t, alphabet_type>>>;

    //!\brief Befriend the base type.
    friend base_t;

    //!\brief The type of the alphabet character.
    using char_type = detail::valid_template_spec_or_t<char, alphabet_char_t, alphabet_type>;

    //!\brief The type of the Phred score.
    using phred_type = detail::valid_template_spec_or_t<int8_t, alphabet_phred_t, alphabet_type>;

private:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_proxy() noexcept = default;                          //!< Defaulted.
    constexpr alphabet_proxy(alphabet_proxy const &) = default;             //!< Defaulted.
    constexpr alphabet_proxy(alphabet_proxy &&) = default;                  //!< Defaulted.
    constexpr alphabet_proxy & operator=(alphabet_proxy const &) = default; //!< Defaulted.
    constexpr alphabet_proxy & operator=(alphabet_proxy &&) = default;      //!< Defaulted.
    ~alphabet_proxy() = default;                                            //!< Defaulted.

    //!\brief Construction from the emulated type.
    constexpr alphabet_proxy(alphabet_type const a) noexcept
        requires std::is_class_v<alphabet_type>
        : base_t{a}
    {}

    //!\brief Construction from the emulated type.
    constexpr alphabet_proxy(alphabet_type const a) noexcept
        requires (!std::is_class_v<alphabet_type>)
        : base_t{}
    {
        base_t::assign_rank(seqan3::to_rank(a));
    }

    //!\brief Assignment from the emulated type. This function triggers the specialisation in the derived_type.
    constexpr derived_type & operator=(alphabet_type const & c) noexcept
    {
        if constexpr (std::is_class_v<alphabet_type>)
            seqan3::assign_rank_to(seqan3::to_rank(c), static_cast<alphabet_type &>(*this));
        else
            base_t::assign_rank(seqan3::to_rank(c));

        static_cast<derived_type &>(*this).on_update(); // <- this invokes the actual proxy behaviour!
        return static_cast<derived_type &>(*this);
    }

    //!\brief Assignment from any type that the emulated type is assignable from.
    template <typename indirect_assignable_type>
    constexpr derived_type & operator=(indirect_assignable_type const & c) noexcept
        requires weakly_assignable_from<alphabet_type, indirect_assignable_type>
    {
        alphabet_type a{};
        a = c;
        return operator=(a);
    }
    //!\}

    //!\brief Befriend the derived type so it can instantiate.
    friend derived_type;

public:
    //!\brief The alphabet size.
    static constexpr auto alphabet_size = seqan3::alphabet_size<alphabet_type>;

    /*!\name Write functions
     * \brief All of these call the emulated type's write functions and then delegate to
     *        the assignment operator which invokes derived behaviour.
     * \{
     */
    //!\brief Assigns a rank.
    constexpr derived_type & assign_rank(alphabet_rank_t<alphabet_type> const r) noexcept
    {
        alphabet_type tmp{};
        assign_rank_to(r, tmp);
        return operator=(tmp);
    }

    //!\brief Assigns a character.
    constexpr derived_type & assign_char(char_type const c) noexcept
        requires writable_alphabet<alphabet_type>
    {
        alphabet_type tmp{};
        assign_char_to(c, tmp);
        return operator=(tmp);
    }

    //!\brief Assigns a Phred score.
    constexpr derived_type & assign_phred(phred_type const c) noexcept
        requires writable_quality_alphabet<alphabet_type>
    {
        alphabet_type tmp{};
        assign_phred_to(c, tmp);
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
        if constexpr (std::is_class_v<alphabet_type>)
            return *this;
        else
            return assign_rank_to(base_t::to_rank(), alphabet_type{});

        /* Instead of static_cast'ing to the alphabet_type which also considers the constructors of the alphabet_type,
         * we explicitly invoke this operator in various places.
         * This prevents errors associated with using alphabet_type's constructors.
         *
         * This is one of error cases:
         * The tuple composite seqan3::qualified returns a component_proxy which inherits from alphabet_proxy_base.
         * The qualified alphabet itself inherits from phred_base.
         * Now when accessing get<1>(seq_qual_alph) we want to call to_phred at some point because we want the quality,
         * therefore the to_phred function from alphabet_proxy is called, but this function did a static_cast to the
         * derived type which is calling the constructor from phred_base. Unfortunately now, the generic phred_base
         * constructor uses `assign_phred_to(to_phred(other), static_cast<derived_type &>(*this))`; (here) which again
         * tries to call to_phred of the alphabet_proxy => infinite loop :boom:
         */
    }

    //!\brief Implicit conversion to types that the emulated type is convertible to.
    template <typename other_t>
        requires (!std::is_class_v<alphabet_type>) && std::convertible_to<alphabet_type, other_t>
    constexpr operator other_t() const noexcept
    {
        return operator alphabet_type();
    }

    //!\brief Returns the rank.
    constexpr auto to_rank() const noexcept
    {
        return seqan3::to_rank(operator alphabet_type());
    }

    //!\brief Returns the character.
    constexpr auto to_char() const noexcept
        requires alphabet<alphabet_type>
    {
        return seqan3::to_char(operator alphabet_type());
    }

    //!\brief Returns the Phred score.
    constexpr auto to_phred() const noexcept
        requires quality_alphabet<alphabet_type>
    {
        return seqan3::to_phred(operator alphabet_type());
    }

    //!\brief Returns the complement.
    constexpr alphabet_type complement() const noexcept
        requires nucleotide_alphabet<alphabet_type>
    {
        return seqan3::complement(operator alphabet_type());
    }

    //!\brief Delegate to the emulated type's validator.
    static constexpr bool char_is_valid(char_type const c) noexcept
        requires writable_alphabet<alphabet_type>
    {
        return char_is_valid_for<alphabet_type>(c);
    }
    //!\}

    /*!\name Comparison operators
     * \brief These are only required if the emulated type allows comparison with types it is not convertible to,
     *        e.g. seqan3::alphabet_variant.
     * \{
     */

private:
    //!\brief work around a gcc bug that disables short-circuiting of operator&& in an enable_if_t of a friend function
    template <typename t>
    static constexpr bool is_alphabet_comparable_with =
        !std::is_same_v<derived_type, t> && detail::weakly_equality_comparable_with<alphabet_type, t>;

public:
    //!\brief Allow (in-)equality comparison with types that the emulated type is comparable with.
    template <typename t>
    friend constexpr auto operator==(derived_type const lhs,
                                     t const rhs) noexcept -> std::enable_if_t<is_alphabet_comparable_with<t>, bool>
    {
        return (lhs.operator alphabet_type() == rhs);
    }

    //!\brief Allow (in-)equality comparison with types that the emulated type is comparable with.
    template <typename t>
    friend constexpr auto
    operator==(t const lhs, derived_type const rhs) noexcept -> std::enable_if_t<is_alphabet_comparable_with<t>, bool>
    {
        return (rhs == lhs);
    }

    //!\brief Allow (in-)equality comparison with types that the emulated type is comparable with.
    template <typename t>
    friend constexpr auto operator!=(derived_type const lhs,
                                     t const rhs) noexcept -> std::enable_if_t<is_alphabet_comparable_with<t>, bool>
    {
        return !(lhs == rhs);
    }

    //!\brief Allow (in-)equality comparison with types that the emulated type is comparable with.
    template <typename t>
    friend constexpr auto
    operator!=(t const lhs, derived_type const rhs) noexcept -> std::enable_if_t<is_alphabet_comparable_with<t>, bool>
    {
        return (rhs != lhs);
    }
    //!\}
};

} // namespace seqan3
