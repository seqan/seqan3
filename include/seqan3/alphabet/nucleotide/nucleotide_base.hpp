// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::nucleotide_base.
 */

#pragma once

#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that refines seqan3::alphabet_base and is used by the nucleotides.
 * \ingroup nucleotide
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 *
 * \details
 *
 * You can use this class to define your own nucleotide alphabet, but types are not required to be based on it to model
 * seqan3::nucleotide_concept, it is purely a way to avoid code duplication.
 *
 * The derived type needs to define only the following table as static member variable:
 *
 *   * `static std::array<THAT_TYPE, value_size> complement_table` that defines for every possible rank value
 *     the corresponding complement.
 */
template <typename derived_type, auto size>
class nucleotide_base : public alphabet_base<derived_type, size, char>
{
private:
    //!\brief Type of the base class.
    using base_t = alphabet_base<derived_type, size, char>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr nucleotide_base() noexcept : base_t{} {}
    constexpr nucleotide_base(nucleotide_base const &) = default;
    constexpr nucleotide_base(nucleotide_base &&) = default;
    constexpr nucleotide_base & operator=(nucleotide_base const &) = default;
    constexpr nucleotide_base & operator=(nucleotide_base &&) = default;
    ~nucleotide_base() = default;
    //!\}

    //! Befriend the derived_type so it can instantiate.
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
    //!\brief Allow explicit construction from any other nucleotide type and convert via the character representation.
    template <typename other_nucl_type>
    //!\cond
        requires !std::Same<nucleotide_base, other_nucl_type> &&
                 !std::Same<derived_type, other_nucl_type> &&
                 nucleotide_concept<other_nucl_type>
    //!\endcond
    explicit constexpr nucleotide_base(other_nucl_type const & other) noexcept
    {
        using seqan3::to_rank;
        static_cast<derived_type &>(*this) =
            detail::convert_through_char_representation<derived_type, other_nucl_type>[to_rank(other)];
    }
    //!\}

    /*!\name Read functions
     * \{
     */

    /*!\brief Return the complement of the letter.
     *
     * \details
     *
     * See \ref nucleotide for the actual values.
     *
     * Satisfies the seqan3::nucleotide_concept::complement() requirement via the seqan3::complement() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr derived_type complement() const noexcept
    {
        return derived_type::complement_table[to_rank()];
    }
    //!\}
};

} // namespace seqan3
