// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::nucleotide_alphabet.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/std/concepts>

// ============================================================================
// complement()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename ...args_t>
void complement(args_t ...) = delete;

//!\brief Functor definition for seqan3::complement.
struct complement_fn
{
public:
    SEQAN3_CPO_IMPL(2, seqan3::custom::alphabet<decltype(v)>::complement(v)) // explicit customisation
    SEQAN3_CPO_IMPL(1, complement(v)                                       ) // ADL
    SEQAN3_CPO_IMPL(0, v.complement()                                      ) // member

public:
    //!\brief Operator definition.
    template <typename nucleotide_t>
    //!\cond
        requires requires (nucleotide_t const nucl)
        {
            { impl(priority_tag<2>{}, nucl) };
            requires noexcept(impl(priority_tag<2>{}, nucl));
            requires std::same_as<nucleotide_t, decltype(impl(priority_tag<2>{}, nucl))>;
        }
    //!\endcond
    constexpr nucleotide_t operator()(nucleotide_t const nucl) const noexcept
    {
        return impl(priority_tag<2>{}, nucl);
    }
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\name Function objects (Nucleotide)
 * \{
 */

/*!\brief Return the complement of a nucleotide object.
 * \tparam your_type Type of the argument.
 * \param  nucl      The nucleotide object for which you want to receive the complement.
 * \returns The complement character of `nucl`, e.g. 'C' for 'G'.
 * \ingroup nucleotide
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A static member function `complement(your_type const a)` of the class `seqan3::custom::alphabet<your_type>`.
 *   2. A free function `complement(your_type const a)` in the namespace of your type (or as `friend`).
 *   3. A member function called `complement()`.
 *
 * Functions are only considered for one of the above cases if they are marked `noexcept` (`constexpr` is not required,
 * but recommended) and if the returned type is `your_type`.
 *
 * Every nucleotide alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/nucleotide/complement_fn.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto complement = detail::adl_only::complement_fn{};
//!\}

// ============================================================================
// nucleotide_alphabet concept
// ============================================================================

/*!\interface seqan3::nucleotide_alphabet <>
 * \extends seqan3::alphabet
 * \brief A concept that indicates whether an alphabet represents nucleotides.
 * \ingroup nucleotide
 *
 * \details
 *
 * In addition to the requirements for seqan3::alphabet, the nucleotide_alphabet introduces
 * a requirement for a complement function: seqan3::nucleotide_alphabet::complement.
 *
 * ### Requirements
 *
 *   1. `t` shall model seqan3::alphabet
 *   2. seqan3::complement needs to be defined for objects of type `t`
 *
 * See the documentation pages for the respective requirements.
 *
 * ### Related types
 *
 * If a given type `t` models this concept, the following types typically do so, as well:
 *
 *   * `t &`
 *   * `t const`
 *   * `t const &`
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT nucleotide_alphabet = alphabet<t> && requires (t val)
{
    { seqan3::complement(val) };
};
//!\endcond

} // namespace seqan3
