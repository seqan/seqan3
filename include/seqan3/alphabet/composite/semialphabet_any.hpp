// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::semialphabet_any.
 */

#pragma once

#include <seqan3/alphabet/alphabet_base.hpp>

namespace seqan3
{

/*!\brief A semi-alphabet that type erases all other semi-alphabets of the same size.
 * \ingroup alphabet_composite
 * \implements seqan3::semialphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 *
 * \details
 * This alphabet provides a generic representation for different alphabets of the same size by erasing
 * the type of the alphabet it is constructed from. This enables the usage of a single type although
 * assigning different alphabets. The semialphabet_any can also be converted to any other (semi-)alphabet
 * of the same size.
 * It is therefore possible to convert the semialphabet_any into an alphabet type that is not the original
 * alphabet type. However, this should either be avoided or used with care as no warnings are given when attempting
 * to convert the semialphabet_any into a type that is not comparable to the original alphabet type.
 * The main advantage of using this data structure is to reduce instantiation of templates when using multiple
 * alphabets of the same size and either their character representation is not important or they are reified at
 * a later point in the program.
 *
 * \see https://en.wikipedia.org/wiki/Type_erasure
 * \see https://en.wikipedia.org/wiki/Reification_(computer_science)
 *
 * ### Example
 * \include test/snippet/alphabet/composite/semialphabet_any.cpp
 *
 * \stableapi{Since version 3.1.}
 */
template <size_t size>
class semialphabet_any : public alphabet_base<semialphabet_any<size>, size, void>
{
private:
    //!\brief Type of the base class.
    using base_t = alphabet_base<semialphabet_any<size>, size, void>;

    //!\brief Befriend the base class.
    friend base_t;

public:
    using base_t::assign_rank;
    using base_t::to_rank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr semialphabet_any() noexcept = default;                                     //!< Defaulted.
    constexpr semialphabet_any(semialphabet_any const &) noexcept = default;             //!< Defaulted.
    constexpr semialphabet_any(semialphabet_any &&) noexcept = default;                  //!< Defaulted.
    constexpr semialphabet_any & operator=(semialphabet_any const &) noexcept = default; //!< Defaulted.
    constexpr semialphabet_any & operator=(semialphabet_any &&) noexcept = default;      //!< Defaulted.
    ~semialphabet_any() noexcept = default;                                              //!< Defaulted.

    /*!\brief Construct semialphabet_any from alphabet of the same size
     * \details
     * \stableapi{Since version 3.1.}
     */
    template <typename other_alph_t>
        requires (!std::same_as<other_alph_t, semialphabet_any>)
              && semialphabet<other_alph_t> && (alphabet_size<other_alph_t> == size)
    explicit semialphabet_any(other_alph_t const other)
    {
        assign_rank(seqan3::to_rank(other));
    }
    //!\}

    /*!\brief Enable conversion of semialphabet_any into other (semi-)alphabet of the same size
     * \details
     * \stableapi{Since version 3.1.}
     */
    template <typename other_alph_t>
        requires (!std::same_as<other_alph_t, semialphabet_any>)
              && semialphabet<other_alph_t> && (alphabet_size<other_alph_t> == size) && std::regular<other_alph_t>
    explicit operator other_alph_t() const
    {
        other_alph_t other{};
        assign_rank_to(to_rank(), other);
        return other;
    }
};

} // namespace seqan3
