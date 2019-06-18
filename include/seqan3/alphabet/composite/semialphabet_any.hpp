// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains seqan3::semialphabet_any.
 */

#pragma once

#include <seqan3/alphabet/detail/alphabet_base.hpp>

namespace seqan3
{

/*!\brief A semi-alphabet that type erases all other semi-alphabets of the same size.
 * \ingroup composite
 * \implements seqan3::Semialphabet
 * \implements seqan3::TriviallyCopyable
 * \implements seqan3::StandardLayout
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
 * \snippet test/snippet/alphabet/composite/semialphabet_any.cpp example
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
    using base_t::to_rank;
    using base_t::assign_rank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr semialphabet_any()                                     noexcept = default; //!< Defaulted.
    constexpr semialphabet_any(semialphabet_any const &)             noexcept = default; //!< Defaulted.
    constexpr semialphabet_any(semialphabet_any &&)                  noexcept = default; //!< Defaulted.
    constexpr semialphabet_any & operator=(semialphabet_any const &) noexcept = default; //!< Defaulted.
    constexpr semialphabet_any & operator=(semialphabet_any &&)      noexcept = default; //!< Defaulted.
    ~semialphabet_any()                                              noexcept = default; //!< Defaulted.

    //!\brief Construct semialphabet_any from alphabet of the same size
    template <Semialphabet other_alph_t>
    //!\cond
        requires (alphabet_size<other_alph_t> == size)
    //!\endcond
    explicit semialphabet_any(other_alph_t const other)
    {
        assign_rank(seqan3::to_rank(other));
    }
    //!\}

    //!\brief Enable conversion of semialphabet_any into other (semi-)alphabet of the same size
    template <Semialphabet other_alph_t>
    //!\cond
        requires (alphabet_size<other_alph_t> == size)
    //!\endcond
    explicit operator other_alph_t() const
    {
        other_alph_t other{};
        assign_rank_to(to_rank(), other);
        return other;
    }
};

} // namespace seqan3
