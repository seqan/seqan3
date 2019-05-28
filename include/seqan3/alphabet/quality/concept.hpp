// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Quality alphabet concept.
 */

#pragma once

#include <iostream>
#include <string>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{

// ============================================================================
// Member exposure for the seqan3::QualityAlphabet
// ============================================================================

/*!\name Helpers for seqan3::QualityAlphabet
 * \brief These functions and metafunctions expose member variables and types so
 * that the type can model the seqan3::QualityAlphabet.
 * \ingroup quality
 * \{
 */
template <typename alphabet_type>
struct underlying_phred
{};

/*!\brief The internal phred type.
 * \ingroup quality
 * \tparam alphabet_type The type of alphabet. Must model the seqan3::QualityAlphabet.
 *
 * The underlying_phred type requires the QualityAlphabet.
 */
template <typename alphabet_with_member_type>
//!\cond
    requires requires () { typename alphabet_with_member_type::phred_type; }
//!\endcond
struct underlying_phred<alphabet_with_member_type>
{
    //!\brief The underlying phred data type.
    using type = typename alphabet_with_member_type::phred_type;
};

/*!\brief The internal phred type.
 * \ingroup quality
 * \tparam alphabet_type The type of alphabet. Must model the seqan3::QualityAlphabet.
 *
 * The underlying_phred type requires the QualityAlphabet.
 */
template <typename alphabet_type>
using underlying_phred_t = typename underlying_phred<alphabet_type>::type;

} // namespace seqan3

// ============================================================================
// assign_phred_to()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::assign_phred_to.
//!\ingroup quality
struct assign_phred_to_fn
{
private:
    SEQAN3_CPO_IMPL(1, (assign_phred_to(args..., v) )) // ADL
    SEQAN3_CPO_IMPL(0, (v.assign_phred(args...)     )) // member

public:
    //!\brief Operator definition for lvalues.
    template <typename alph_t>
    //!\cond
        requires requires (seqan3::underlying_phred_t<alph_t> const p, alph_t & a)
            { { impl(priority_tag<1>{}, a, p) }; }
    //!\endcond
    constexpr alph_t & operator()(seqan3::underlying_phred_t<alph_t> const p, alph_t & a) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<1>{}, a, p)),
            "Only overloads that are marked noexcept are picked up by seqan3::assign_phred_to().");
        static_assert(std::Same<alph_t &, decltype(impl(priority_tag<1>{}, a, p))>,
            "The return type of your assign_phred_to() implementation must be 'alph_t &'.");

        return impl(priority_tag<1>{}, a, p);
    }

    //!\brief Operator definition for rvalues.
    template <typename alph_t>
    //!\cond
        requires requires (seqan3::underlying_phred_t<alph_t> const p, alph_t & a)
            { { impl(priority_tag<1>{}, a, p) }; } && (!std::is_lvalue_reference_v<alph_t>)
    //!\endcond
    constexpr alph_t operator()(seqan3::underlying_phred_t<alph_t> const p, alph_t && a) const noexcept
    {
        return (*this)(p, a); // call above function but return by value
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Assign a phred score to a quality alphabet object.
 * \tparam your_type The type of the target object. Must model the seqan3::QualityAlphabet.
 * \param  chr       The phred score being assigned; must be of the seqan3::underlying_phred_t of the target object.
 * \returns Reference to `alph` if `alph` was given as lvalue, otherwise a copy.
 * \ingroup quality
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for two possible implementations (in this order):
 *
 *   1. A free function `assign_phred_to(phred_type const chr, your_type & a)` in the namespace of your type
 *      (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be `your_type &`.
 *   2. A member function called `assign_phred(phred_type const chr)` (not `assign_phred_to`).
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `your_type &`.
 *
 * Every writable quality alphabet type must provide one of the above. *Note* that temporaries of `your_type`
 * are handled by this function object and **do not** require an additional overload.
 *
 * ### Customisation point
 *
 * This is a customisation point. To specify the behaviour for your own alphabet type,
 * simply provide one of the two functions specified above.
 */
inline constexpr auto assign_phred_to = detail::adl::only::assign_phred_to_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// to_phred()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::to_phred.
struct to_phred_fn
{
private:
    SEQAN3_CPO_IMPL(1, to_phred(v)  ) // ADL
    SEQAN3_CPO_IMPL(0, v.to_phred() ) // member

public:
    //!\brief Operator definition.
    template <typename alph_t>
    //!\cond
        requires requires (alph_t const chr) { { impl(priority_tag<1>{}, chr) }; }
    //!\endcond
    constexpr auto operator()(alph_t const chr) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<1>{}, chr)),
            "Only overloads that are marked noexcept are picked up by seqan3::to_phred().");
        static_assert(std::Constructible<size_t, decltype(impl(priority_tag<1>{}, chr))>,
            "The return type of your to_phred() implementation must be convertible to size_t.");

        return impl(priority_tag<1>{}, chr);
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief The public getter function for the phred representation of a quality score.
 * \tparam your_type The type of alphabet. Must model the seqan3::QualityAlphabet.
 * \param  chr       The quality value to convert into the phred score.
 * \returns the phred representation of a quality score.
 * \ingroup quality
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for two possible implementations (in this order):
 *
 *   1. A free function `to_phred(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be of the respective phred representation (usually a small integral type).
 *   2. A member function called `to_phred()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be of
 *      the respective rank representation.
 *
 * Every quality alphabet type must provide one of the above.
 *
 * ### Customisation point
 *
 * This is a customisation point. To specify the behaviour for your own alphabet type,
 * simply provide one of the two functions specified above.
 */
inline constexpr auto to_phred = detail::adl::only::to_phred_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// seqan3::QualityAlphabet
// ============================================================================

namespace seqan3
{

/*!\interface seqan3::QualityAlphabet <>
 * \extends seqan3::Alphabet
 * \brief A concept that indicates whether an alphabet represents quality scores.
 * \ingroup quality
 *
 * In addition to the requirements for seqan3::Alphabet, the
 * QualityAlphabet introduces a requirement for conversion functions from and to
 * a Phred score.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and
 * metafunctions. Types that satisfy this concept are shown as "implementing
 * this interface".
 */
//!\cond
template<typename q>
SEQAN3_CONCEPT QualityAlphabet = requires(q quality)
{
    requires Alphabet<q>;

    { to_phred(quality) } -> const typename q::phred_type;
    typename underlying_phred<q>::type;
};
//!\endcond

// ============================================================================
// seqan3::WritableQualityAlphabet
// ============================================================================

/*!\interface seqan3::QualityAlphabet <>
 * \extends seqan3::Alphabet
 * \brief A concept that indicates whether a writable alphabet represents quality scores.
 * \ingroup quality
 *
 * In addition to the requirements for seqan3::WritableAlphabet, the seqan3::WritableQualityAlphabet
 * introduces the requirements of seqan3::QualityAlphabet.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and
 * metafunctions. Types that satisfy this concept are shown as "implementing
 * this interface".
 */
//!\cond
template<typename q>
SEQAN3_CONCEPT WritableQualityAlphabet = requires(q quality)
{
    requires WritableAlphabet<q>;
    requires QualityAlphabet<q>;

    { assign_phred_to(typename q::rank_type{}, quality) } -> q;
};
//!\endcond

} // namespace seqan3
