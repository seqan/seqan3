// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::nucleotide_alphabet.
 */

#pragma once

#include <concepts>

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// complement()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void complement(args_t...) = delete;

//!\brief seqan3::detail::customisation_point_object (CPO) definition for seqan3::complement.
//!\ingroup alphabet_nucleotide
struct complement_cpo : public detail::customisation_point_object<complement_cpo, 2>
{
    //!\brief CRTP base class seqan3::detail::customisation_point_object.
    using base_t = detail::customisation_point_object<complement_cpo, 2>;
    //!\brief Only this class is allowed to import the constructors from #base_t. (CRTP safety idiom)
    using base_t::base_t;

    /*!\brief CPO overload (check 1 out of 3): explicit customisation via `seqan3::custom::alphabet`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet the complement is returned from.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<2>, alphabet_t && alphabet)(
        /*return*/ seqan3::custom::alphabet<alphabet_t>::complement(std::forward<alphabet_t>(alphabet)) /*;*/
    );

    /*!\brief CPO overload (check 2 out of 3): argument dependent lookup (ADL), i.e. `complement(alphabet)`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet the complement is returned from.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<1>, alphabet_t && alphabet)(
        /*return*/ complement(std::forward<alphabet_t>(alphabet)) /*;*/
    );

    /*!\brief CPO overload (check 3 out of 3): member access, i.e. `alphabet.complement()`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet the complement is returned from.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<0>, alphabet_t && alphabet)(
        /*return*/ std::forward<alphabet_t>(alphabet).complement() /*;*/
    );
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
 * \ingroup alphabet_nucleotide
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
 * \include test/snippet/alphabet/nucleotide/complement_cpo.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 *
 * \experimentalapi{Implementation 2 (free function) is not stable.}
 *
 * \stableapi{Since version 3.1. The name seqan3::complement, Implementation 1,
 *            and Implementation 3 are stable and will not change.}
 */
inline constexpr auto complement = detail::adl_only::complement_cpo{};
//!\}

// ============================================================================
// nucleotide_alphabet concept
// ============================================================================

/*!\interface seqan3::nucleotide_alphabet <>
 * \extends seqan3::alphabet
 * \brief A concept that indicates whether an alphabet represents nucleotides.
 * \ingroup alphabet_nucleotide
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
 *
 * \stableapi{Since version 3.1.}
 */
//!\cond
template <typename t>
concept nucleotide_alphabet = alphabet<t> && requires (t val) {
    { seqan3::complement(val) };
};
//!\endcond

} // namespace seqan3
