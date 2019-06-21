// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Extends a given alphabet with the mask alphabet.
 */

#pragma once

#include <seqan3/alphabet/mask/all.hpp>
#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/core/char_operations/transform.hpp>

namespace seqan3
{
/*!\brief Implementation of a masked composite, which extends a given alphabet
 * with a mask.
 * \ingroup mask
 * \implements seqan3::WritableAlphabet
 * \if DEV \implements seqan3::detail::WritableConstexprAlphabet \endif
 * \implements seqan3::TriviallyCopyable
 * \implements seqan3::StandardLayout
 * \implements std::Regular
 *
 * \tparam sequence_alphabet_t Type of the first letter; must satisfy seqan3::Semialphabet.
 * \tparam mask_t Types of masked letter; must satisfy seqan3::Semialphabet, defaults to seqan3::mask.
 *
 * \details
 * The masked composite represents a seqan3::alphabet_tuple_base of any given alphabet with the
 * masked alphabet. It allows one to specify which portions of a sequence should be masked,
 * without losing additional information by replacing the sequence directly.
 *
 * \include test/snippet/alphabet/mask/masked.cpp
 */
 template <typename sequence_alphabet_t>
//!\cond
    requires WritableAlphabet<sequence_alphabet_t>
//!\endcond
class masked : public alphabet_tuple_base<masked<sequence_alphabet_t>, sequence_alphabet_t, mask>
{
private:
    //!\brief The base type.
    using base_type = alphabet_tuple_base<masked<sequence_alphabet_t>, sequence_alphabet_t, mask>;

public:
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = alphabet_char_t<sequence_alphabet_type>;

    using base_type::alphabet_size;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr masked() : base_type{} {}                      //!< Defaulted.
    constexpr masked(masked const &) = default;              //!< Defaulted.
    constexpr masked(masked &&) = default;                   //!< Defaulted.
    constexpr masked & operator =(masked const &) = default; //!< Defaulted.
    constexpr masked & operator =(masked &&) = default;      //!< Defaulted.
    ~masked() = default;                                     //!< Defaulted.

    using base_type::base_type; // Inherit non-default constructors

    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked(component_type const alph) {} ))
    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked(indirect_component_type const alph) {} ))
    //!\copydoc alphabet_tuple_base::operator=(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked & operator=(component_type const alph) {} ))
    //!\copydoc alphabet_tuple_base::operator=(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr masked & operator=(indirect_component_type const alph) {} ))
    //!\}

    // Inherit operators from base
    using base_type::operator=;
    using base_type::operator==;
    using base_type::operator!=;
    using base_type::operator>=;
    using base_type::operator<=;
    using base_type::operator<;
    using base_type::operator>;

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    constexpr masked & assign_char(char_type const c) noexcept
    {
        seqan3::assign_char_to(c, get<0>(*this));
        seqan3::assign_rank_to(is_lower(c), get<1>(*this));
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return a character.
    constexpr char_type to_char() const noexcept
    {
        if (seqan3::to_rank(get<1>(*this)))
        {
            return to_lower(seqan3::to_char(get<0>(*this)));
        }
        else
        {
            return seqan3::to_char(get<0>(*this));
        }
    }
    //!\}

    /*!\brief Validate whether a character value has a one-to-one mapping to an alphabet value.
     *
     * \details
     *
     * Satisfies the seqan3::Semialphabet::char_is_valid_for() requirement via the seqan3::char_is_valid_for()
     * wrapper.
     *
     * Default implementation: True for all character values that are reproduced by #to_char() after being assigned
     * to the alphabet.
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
        return masked{}.assign_char(c).to_char() == c;
    }
};

//!\brief Type deduction guide enables usage of masked without specifying template args.
//!\relates masked
template <typename sequence_alphabet_type>
masked(sequence_alphabet_type &&, mask const &)
    -> masked<std::decay_t<sequence_alphabet_type>>;
} //namespace seqan3
