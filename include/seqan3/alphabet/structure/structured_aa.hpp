// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides the composite of aminoacid with structure alphabets.
 */

#pragma once

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/structure/dssp9.hpp>

namespace seqan3
{

/*!\brief A seqan3::alphabet_tuple_base that joins an aminoacid alphabet with a protein structure alphabet.
 * \ingroup structure
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \tparam sequence_alphabet_t Type of the first aminoacid letter; must model seqan3::writable_alphabet,
 * seqan3::aminoacid_alphabet and satisfy the requirements on arguments by seqan3::alphabet_tuple_base.
 * \tparam structure_alphabet_t Types of further structure letters; must model seqan3::writable_alphabet and satisfy
 * the requirements on arguments by seqan3::alphabet_tuple_base.
 *
 * This composite pairs an aminoacid alphabet with a structure alphabet. The rank values
 * correpsond to numeric values in the size of the composite, while the character values
 * are taken from the sequence alphabet and the structure annotation is taken from the structure
 * alphabet.
 *
 * As with all `seqan3::alphabet_tuple_base` s you may access the individual alphabet letters in
 * regular c++ tuple notation, i.e. `get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * \include test/snippet/alphabet/structure/structured_aa.cpp
 *
 * This seqan3::alphabet_tuple_base itself fulfills seqan3::alphabet.
 */
template <writable_alphabet sequence_alphabet_t = aa27, writable_alphabet structure_alphabet_t = dssp9>
//!\cond
    requires (!std::is_reference_v<sequence_alphabet_t>) && (!std::is_reference_v<structure_alphabet_t>)
//!\endcond
class structured_aa :
    public alphabet_tuple_base<structured_aa<sequence_alphabet_t, structure_alphabet_t>,
                               sequence_alphabet_t, structure_alphabet_t>
{
private:
    //!\brief The base type.
    using base_type = alphabet_tuple_base<structured_aa<sequence_alphabet_t, structure_alphabet_t>,
                                          sequence_alphabet_t, structure_alphabet_t>;
public:
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;
    //!\brief Second template parameter as member type.
    using structure_alphabet_type = structure_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = alphabet_char_t<sequence_alphabet_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr structured_aa()                                    noexcept = default; //!< Defaulted.
    constexpr structured_aa(structured_aa const &)               noexcept = default; //!< Defaulted.
    constexpr structured_aa(structured_aa &&)                    noexcept = default; //!< Defaulted.
    constexpr structured_aa & operator =(structured_aa const &)  noexcept = default; //!< Defaulted.
    constexpr structured_aa & operator =(structured_aa &&)       noexcept = default; //!< Defaulted.
    ~structured_aa()                                             noexcept = default; //!< Defaulted.

    using base_type::base_type; // Inherit non-default constructors

#if SEQAN3_DOXYGEN_ONLY(1)0
    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(component_type const alph)
    template <typename component_type>
    constexpr structured_aa(component_type const alph) {}

    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(indirect_component_type const alph)
    template <typename indirect_component_type>
    constexpr structured_aa(indirect_component_type const alph) {}

    //!\copydoc alphabet_tuple_base::operator=(component_type const alph)
    template <typename component_type>
    constexpr structured_aa & operator=(component_type const alph) {}

    //!\copydoc alphabet_tuple_base::operator=(indirect_component_type const alph)
    template <typename indirect_component_type>
    constexpr structured_aa & operator=(indirect_component_type const alph) {}
#endif

    //!\brief Inherit operators from base
    using base_type::operator=;
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a nucleotide character. This modifies the internal sequence letter.
    constexpr structured_aa & assign_char(char_type const c) noexcept
    {
        seqan3::assign_char_to(c, get<0>(*this));
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return a character. This reads the internal sequence letter.
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }
    //!\}

    //!\brief Validate whether a character is valid in the sequence alphabet.
    static constexpr bool char_is_valid(char_type const c) noexcept
    {
        return char_is_valid_for<sequence_alphabet_type>(c);
    }
};

//!\brief Type deduction guide enables usage of structured_aa without specifying template args.
//!\relates structured_aa
template <typename sequence_alphabet_type, typename structure_alphabet_type>
structured_aa(sequence_alphabet_type &&, structure_alphabet_type &&)
    -> structured_aa<std::decay_t<sequence_alphabet_type>, std::decay_t<structure_alphabet_type>>;

} // namespace seqan3
