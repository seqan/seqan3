// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::NucleotideAlphabet.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/std/concepts>

// ============================================================================
// complement()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::complement.
struct complement_fn
{
private:
    SEQAN3_CPO_IMPL(1, complement(v)                       )    // ADL
    SEQAN3_CPO_IMPL(0, v.complement()                      )    // member

public:
    //!\brief Operator definition.
    template <typename nucleotide_t>
    //!\cond
        requires requires (nucleotide_t const nucl) { { impl(priority_tag<1>{}, nucl) }; }
    //!\endcond
    constexpr auto operator()(nucleotide_t const nucl) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<1>{}, nucl)),
            "Only overloads that are marked noexcept are picked up by seqan3::complement().");
        static_assert(std::Same<nucleotide_t, decltype(impl(priority_tag<1>{}, nucl))>,
            "The return type of your complement() implementation must be 'nucleotide_t'.");

        return impl(priority_tag<1>{}, nucl);
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Return the complement of a nucleotide object.
 * \tparam your_type Type of the argument.
 * \param  nucl      The nucleotide object.
 * \returns The complement character of `nucl`, e.g. 'C' for 'G'.
 * \ingroup nucleotide
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for two possible implementations (in this order):
 *
 *   1. A free function `complement(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be `your_type`.
 *   2. A member function called `complement()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `your_type`.
 *
 * Every nucleotide alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/nucleotide/complement_fn.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point. To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto complement = detail::adl::only::complement_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// NucleotideAlphabet concept
// ============================================================================

namespace seqan3
{

/*!\interface seqan3::NucleotideAlphabet <>
 * \extends seqan3::Alphabet
 * \brief A concept that indicates whether an alphabet represents nucleotides.
 * \ingroup nucleotide
 *
 * In addition to the requirements for seqan3::Alphabet, the NucleotideAlphabet introduces
 * a requirement for a complement function: seqan3::NucleotideAlphabet::complement.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT NucleotideAlphabet = Alphabet<type> && requires (type v, remove_cvref_t<type> c)
{
    requires std::Same<decltype(complement(v)), decltype(c)>;
};
//!\endcond

/*!\name Requirements for seqan3::NucleotideAlphabet
 * \brief You can expect these functions on all types that implement seqan3::NucleotideAlphabet.
 * \{
 */
/*!\fn nucleotide_type seqan3::complement(nucleotide_type const alph)
 * \brief Returns the alphabet letter's complement value.
 * \relates seqan3::NucleotideAlphabet
 * \param alph The alphabet letter for whom you wish to receive the complement.
 * \returns The letter's complement, e.g. 'T' for 'A'.
 * \details
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
} // namespace seqan3
