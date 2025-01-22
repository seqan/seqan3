// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides the composite of nucleotide with structure alphabets.
 */

#pragma once

#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/structure/concept.hpp>

namespace seqan3
{

/*!\brief A seqan3::alphabet_tuple_base that joins a nucleotide alphabet with an RNA structure alphabet.
 * \ingroup alphabet_structure
 * \implements seqan3::rna_structure_alphabet
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \tparam sequence_alphabet_t Type of the first letter; must model seqan3::writable_alphabet,
 * seqan3::nucleotide_alphabet and satisfy the requirements on arguments by seqan3::alphabet_tuple_base.
 * \tparam structure_alphabet_t Types of further letters; must model seqan3::writable_alphabet,
 * seqan3::rna_structure_alphabet and satisfy the requirements on arguments by seqan3::alphabet_tuple_base.
 *
 * This composite pairs a nucleotide alphabet with a structure alphabet. The rank values
 * correspond to numeric values in the size of the composite, while the character values
 * are taken from the sequence alphabet and the structure annotation is taken from the structure
 * alphabet.
 *
 * As with all `seqan3::alphabet_tuple_base` s you may access the individual alphabet letters in
 * regular c++ tuple notation, i.e. `get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * \include test/snippet/alphabet/structure/structured_rna.cpp
 *
 * This seqan3::alphabet_tuple_base itself models both seqan3::nucleotide_alphabet and seqan3::rna_structure_alphabet.
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
template <nucleotide_alphabet sequence_alphabet_t, rna_structure_alphabet structure_alphabet_t>
    requires writable_alphabet<sequence_alphabet_t> && writable_alphabet<structure_alphabet_t>
class structured_rna :
    public alphabet_tuple_base<structured_rna<sequence_alphabet_t, structure_alphabet_t>,
                               sequence_alphabet_t,
                               structure_alphabet_t>
{
private:
    //!\brief The base type.
    using base_type = alphabet_tuple_base<structured_rna<sequence_alphabet_t, structure_alphabet_t>,
                                          sequence_alphabet_t,
                                          structure_alphabet_t>;

public:
    /*!\brief First template parameter as member type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using sequence_alphabet_type = sequence_alphabet_t;
    /*!\brief Second template parameter as member type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using structure_alphabet_type = structure_alphabet_t;

    /*!\brief Equals the char_type of sequence_alphabet_type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using char_type = alphabet_char_t<sequence_alphabet_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr structured_rna() noexcept = default;                                   //!< Defaulted.
    constexpr structured_rna(structured_rna const &) noexcept = default;             //!< Defaulted.
    constexpr structured_rna(structured_rna &&) noexcept = default;                  //!< Defaulted.
    constexpr structured_rna & operator=(structured_rna const &) noexcept = default; //!< Defaulted.
    constexpr structured_rna & operator=(structured_rna &&) noexcept = default;      //!< Defaulted.
    ~structured_rna() noexcept = default;                                            //!< Defaulted.

    using base_type::base_type; // Inherit non-default constructors

#if SEQAN3_DOXYGEN_ONLY(1) 0
    /*!\copybrief seqan3::alphabet_tuple_base::alphabet_tuple_base(component_type const alph)
     * \details
     * \sa seqan3::alphabet_tuple_base::alphabet_tuple_base(component_type const alph)
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr structured_rna(component_type const alph) noexcept
    {}

    /*!\copybrief seqan3::alphabet_tuple_base::alphabet_tuple_base(indirect_component_type const alph)
     * \details
     * \sa seqan3::alphabet_tuple_base::alphabet_tuple_base(indirect_component_type const alph)
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr structured_rna(indirect_component_type const alph) noexcept
    {}

    /*!\copybrief seqan3::alphabet_tuple_base::operator=(component_type const alph)
     * \details
     * \sa seqan3::alphabet_tuple_base::operator=(component_type const alph)
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr structured_rna & operator=(component_type const alph) noexcept
    {}

    /*!\copybrief seqan3::alphabet_tuple_base::operator=(indirect_component_type const alph)
     * \details
     * \sa seqan3::alphabet_tuple_base::operator=(indirect_component_type const alph)
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr structured_rna & operator=(indirect_component_type const alph) noexcept
    {}
#endif // SEQAN3_DOXYGEN_ONLY
    //!\}

    // Inherit operators from base
    using base_type::operator=;

    //!\name Write functions
    //!\{

    /*!\brief Assign from a nucleotide character. This modifies the internal sequence letter.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr structured_rna & assign_char(char_type const c) noexcept
    {
        seqan3::assign_char_to(c, get<0>(*this));
        return *this;
    }
    //!\}

    //!\name Read functions
    //!\{

    /*!\brief Return a character. This reads the internal sequence letter.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }

    /*!\brief Return a structured_rna where the sequence letter is converted to its complement.
     * \details
     * See \ref alphabet_nucleotide for the actual values.
     * Satisfies the seqan3::nucleotide_alphabet::complement() requirement via the seqan3::complement() wrapper.
     * The structure letter is not modified.
     * ### Complexity
     * Constant.
     * ### Exceptions
     * Guaranteed not to throw.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr structured_rna complement() const noexcept
    {
        return structured_rna{get<0>(*this).complement(), get<1>(*this)};
    }
    //!\}

    /*!\brief Validate whether a character is valid in the sequence alphabet.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    static constexpr bool char_is_valid(char_type const c) noexcept
    {
        return char_is_valid_for<sequence_alphabet_type>(c);
    }

    //!\name RNA structure properties
    //!\{

    /*!\brief Check whether the character represents a rightward interaction in an RNA structure.
     * \returns True if the letter represents a rightward interaction, False otherwise.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_pair_open() const noexcept
    {
        return get<1>(*this).is_pair_open();
    };

    /*!\brief Check whether the character represents a leftward interaction in an RNA structure.
     * \returns True if the letter represents a leftward interaction, False otherwise.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_pair_close() const noexcept
    {
        return get<1>(*this).is_pair_close();
    };

    /*!\brief Check whether the character represents an unpaired position in an RNA structure.
     * \returns True if the letter represents an unpaired site, False otherwise.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool is_unpaired() const noexcept
    {
        return get<1>(*this).is_unpaired();
    };

    //!\brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions.
    static constexpr uint8_t max_pseudoknot_depth{structure_alphabet_t::max_pseudoknot_depth};

    /*!\brief Get an identifier for a pseudoknotted interaction.
     * \returns The pseudoknot id, if alph denotes an interaction, and no value otherwise.
     * \details
     * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr std::optional<uint8_t> pseudoknot_id() const noexcept
    {
        return get<1>(*this).pseudoknot_id();
    };
    //!\}
};

//!\brief Type deduction guide enables usage of structured_rna without specifying template args.
//!\relates structured_rna
template <typename sequence_alphabet_type, typename structure_alphabet_type>
structured_rna(sequence_alphabet_type &&, structure_alphabet_type &&)
    -> structured_rna<std::decay_t<sequence_alphabet_type>, std::decay_t<structure_alphabet_type>>;

} // namespace seqan3
