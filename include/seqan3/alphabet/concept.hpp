// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Core alphabet concept and free function/metafunction wrappers.
 */

#pragma once

#include <iostream>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/detail/member_exposure.hpp>
#include <seqan3/alphabet/exception.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/customisation_point.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/function.hpp>
#include <seqan3/std/concepts>

// ============================================================================
// to_rank()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::to_rank.
struct to_rank_fn
{
private:
    SEQAN3_CPO_IMPL(2, to_rank(v)                       )    // ADL
    SEQAN3_CPO_IMPL(1, seqan3::adaptation::to_rank(v)   )    // customisation namespace
    SEQAN3_CPO_IMPL(0, v.to_rank()                      )    // member

public:
    //!\brief Operator definition.
    template <typename alph_t>
    //!\cond
        requires requires (alph_t const a) { { impl(priority_tag<2>{}, a) }; }
    //!\endcond
    constexpr auto operator()(alph_t const a) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, a)),
            "Only overloads that are marked noexcept are picked up by seqan3::to_rank().");
        static_assert(std::Constructible<size_t, decltype(impl(priority_tag<2>{}, a))>,
            "The return type of your to_rank() implementation must be convertible to size_t.");

        return impl(priority_tag<2>{}, a);
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Return the rank representation of a (semi-)alphabet object.
 * \tparam your_type Type of the argument.
 * \param  alph      The (semi-)alphabet object.
 * \returns The rank representation; an integral type.
 * \ingroup alphabet
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A free function `to_rank(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be of the respective rank representation (usually a small integral type).
 *   2. A free function `to_rank(your_type const a)` in `namespace seqan3::adaptation`.
 *      Used only by the SeqAn library internally. The same restrictions apply as above.
 *   3. A member function called `to_rank()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      of the respective rank representation.
 *
 * Every (semi-)alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/to_rank.cpp
 *
 * For an example of a full alphabet definition with free function implementations (solution 1. above),
 * see seqan3::Alphabet.
 *
 * ### Customisation point
 *
 * This is a customisation point (TODO link to manual). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto to_rank = detail::adl::only::to_rank_fn{};
//!\}

//!\brief The `rank_type` of the semi-alphabet; defined as the return type of seqan3::to_rank.
//!\ingroup alphabet
template <typename semi_alphabet_type>
//!\cond
    requires requires { { seqan3::to_rank(std::declval<semi_alphabet_type>()) }; }
//!\endcond
using alphabet_rank_t = decltype(seqan3::to_rank(std::declval<semi_alphabet_type>()));

} // namespace seqan3

// ============================================================================
// assign_rank_to()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::assign_rank_to.
//!\ingroup alphabet
struct assign_rank_to_fn
{
private:
    SEQAN3_CPO_IMPL(2, (assign_rank_to(args..., v)                           ))    // ADL
    SEQAN3_CPO_IMPL(1, (seqan3::adaptation::assign_rank_to(args..., v)       ))    // customisation namespace
    SEQAN3_CPO_IMPL(0, (v.assign_rank(args...)                               ))    // member

public:
    //!\brief Operator definition for lvalues.
    template <typename alph_t>
    //!\cond
        requires requires (seqan3::alphabet_rank_t<alph_t> const r, alph_t & a)
            { { impl(priority_tag<2>{}, a, r) }; }
    //!\endcond
    constexpr alph_t & operator()(seqan3::alphabet_rank_t<alph_t> const r, alph_t & a) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, a, r)),
            "Only overloads that are marked noexcept are picked up by seqan3::assign_rank_to().");
        static_assert(std::Same<alph_t &, decltype(impl(priority_tag<2>{}, a, r))>,
            "The return type of your assign_rank_to() implementation must be 'alph_t &'.");

        return impl(priority_tag<2>{}, a, r);
    }

    //!\brief Operator definition for rvalues.
    template <typename alph_t>
    //!\cond
        requires requires (seqan3::alphabet_rank_t<alph_t> const r, alph_t & a)
            { { impl(priority_tag<2>{}, a, r) }; } && (!std::is_lvalue_reference_v<alph_t>)
    //!\endcond
    constexpr alph_t operator()(seqan3::alphabet_rank_t<alph_t> const r, alph_t && a) const noexcept
    {
        return (*this)(r, a); // call above function but return by value
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Assign a rank to an alphabet object.
  *\tparam your_type Type of the target object.
 * \param chr  The rank being assigned; must be of the seqan3::alphabet_rank_t of the target object.
 * \param alph The target object.
 * \returns Reference to `alph` if `alph` was given as lvalue, otherwise a copy.
 * \ingroup alphabet
 *
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A free function `assign_rank_to(rank_type const chr, your_type & a)` in the namespace of your
 *      type (or as `friend`). The function must be marked `noexcept` (`constexpr` is not required,
 *      but recommended) and the return type be `your_type &`.
 *   2. A free function `assign_rank_to(rank_type const chr, your_type & a)` in
 *      `namespace seqan3::adaptation`. Used only by the SeqAn library internally. The same restrictions apply as above.
 *   3. A member function called `assign_rank(rank_type const chr)` (not `assign_rank_to`).
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `your_type &`.
 *
 * Every (semi-)alphabet type must provide one of the above. *Note* that temporaries of `your_type` are handled
 * by this function object and **do not** require an additional overload.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/assign_rank_to.cpp
 *
 * For an example of a full alphabet definition with free function implementations (solution 1. above),
 * see seqan3::Alphabet.
 *
 * ### Customisation point
 *
 * This is a customisation point (TODO link to manual). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto assign_rank_to = detail::adl::only::assign_rank_to_fn{};
//!\}
} // namespace seqan3

// ============================================================================
// to_char()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::to_char.
struct to_char_fn
{
private:
    SEQAN3_CPO_IMPL(2, to_char(v)                       )    // ADL
    SEQAN3_CPO_IMPL(1, seqan3::adaptation::to_char(v)   )    // customisation namespace
    SEQAN3_CPO_IMPL(0, v.to_char()                      )    // member

public:
    //!\brief Operator definition.
    template <typename alph_t>
    //!\cond
        requires requires (alph_t const a) { { impl(priority_tag<2>{}, a) }; }
    //!\endcond
    constexpr decltype(auto) operator()(alph_t const a) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, a)),
            "Only overloads that are marked noexcept are picked up by seqan3::to_char().");

        return impl(priority_tag<2>{}, a);
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Return the char representation of an alphabet object.
 * \tparam your_type Type of the argument.
 * \param  alph      The alphabet object.
 * \returns The char representation; usually `char`.
 * \ingroup alphabet
 *
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A free function `to_char(your_type const a)` in the namespace of your
 *      type (or as `friend`). The function must be marked `noexcept` (`constexpr` is not required,
 *      but recommended) and the return type be of the respective char representation (usually a small integral type).
 *   2. A free function `to_char(your_type const a)` in
 *      `namespace seqan3::adaptation`. Used only by the SeqAn library internally. The same restrictions apply as above.
 *   3. A member function called `to_char()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      of the respective char representation.
 *
 * Every alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/to_char.cpp
 *
 * For an example of a full alphabet definition with free function implementations (solution 1. above),
 * see seqan3::Alphabet.
 *
 * ### Customisation point
 *
 * This is a customisation point (TODO link to manual). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto to_char = detail::adl::only::to_char_fn{};
//!\}

//!\brief The `char_type` of the alphabet; defined as the return type of seqan3::to_char.
//!\ingroup alphabet
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type const a) { { seqan3::to_char(a) }; }
//!\endcond
using alphabet_char_t = decltype(seqan3::to_char(std::declval<alphabet_type const>()));

} // namespace seqan3

// ============================================================================
// assign_char_to()
// ============================================================================

namespace seqan3::adaptation
{
//!\cond
void char_is_valid_for(); // forward
//!\endcond
} // seqan3::adaptation

namespace seqan3::detail::adl::only
{
//!\brief Functor definition for seqan3::assign_char_to.
//!\ingroup alphabet
struct assign_char_to_fn
{
private:
    SEQAN3_CPO_IMPL(2, (assign_char_to(args..., v)                           ))    // ADL
    SEQAN3_CPO_IMPL(1, (seqan3::adaptation::assign_char_to(args..., v)       ))    // customisation namespace
    SEQAN3_CPO_IMPL(0, (v.assign_char(args...)                               ))    // member

public:
    //!\brief Operator definition for lvalues.
    template <typename alph_t>
    //!\cond
        requires requires (seqan3::alphabet_char_t<alph_t> const r, alph_t & a)
            { { impl(priority_tag<2>{}, a, r) }; }
    //!\endcond
    constexpr alph_t & operator()(seqan3::alphabet_char_t<alph_t> const r, alph_t & a) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, a, r)),
            "Only overloads that are marked noexcept are picked up by seqan3::assign_char_to().");
        static_assert(std::Same<alph_t &, decltype(impl(priority_tag<2>{}, a, r))>,
            "The return type of your assign_char_to() implementation must be 'alph_t &'.");

        return impl(priority_tag<2>{}, a, r);
    }

    //!\brief Operator definition for rvalues.
    template <typename alph_t>
    //!\cond
        requires requires (seqan3::alphabet_char_t<alph_t> const r, alph_t & a)
            { { impl(priority_tag<2>{}, a, r) }; } && (!std::is_lvalue_reference_v<alph_t>)
    //!\endcond
    constexpr alph_t operator()(seqan3::alphabet_char_t<alph_t> const r, alph_t && a) const noexcept
    {
        return (*this)(r, a); // call above function but return by value
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Assign a character to an alphabet object.
 * \tparam your_type Type of the target object.
 * \param chr  The character being assigned; must be of the seqan3::alphabet_char_t of the target object.
 * \param alph The target object; its type must model seqan3::Alphabet.
 * \returns Reference to `alph` if `alph` was given as lvalue, otherwise a copy.
 * \ingroup alphabet
 *
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A free function `assign_char_to(char_type const chr, your_type & a)` in the namespace of your
 *      type (or as `friend`). The function must be marked `noexcept` (`constexpr` is not required,
 *      but recommended) and the return type be `your_type &`.
 *   2. A free function `assign_char_to(char_type const chr, your_type & a)` in
 *      `namespace seqan3::adaptation`. Used only by the SeqAn library internally. The same restrictions apply as above.
 *   3. A member function called `assign_char(char_type const chr)` (not `assign_char_to`).
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `your_type &`.
 *
 * Every alphabet type must provide one of the above. *Note* that temporaries of `your_type` are handled
 * by this function object and **do not** require an additional overload.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/assign_char_to.cpp
 *
 * For an example of a full alphabet definition with free function implementations (solution 1. above),
 * see seqan3::Alphabet.
 *
 * ### Customisation point
 *
 * This is a customisation point (TODO link to manual). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto assign_char_to = detail::adl::only::assign_char_to_fn{};
//!\}
} // namespace seqan3

// ============================================================================
// char_is_valid_for()
// ============================================================================

namespace seqan3::adaptation
{
//!\cond
void char_is_valid_for(); // forward
//!\endcond
} // seqan3::adaptation

namespace seqan3::detail::adl::only
{
/*!\brief Functor definition for seqan3::assign_char_strictly_to.
 * \tparam alph_t   The alphabet type being queried.
 * \tparam s_alph_t `alph_t` with cvref removed and possibly wrapped in std::type_identity; never user-provide this!
 * \ingroup alphabet
 */
template <typename alph_t,
          typename s_alph_t = std::conditional_t<std::is_nothrow_default_constructible_v<remove_cvref_t<alph_t>>,
                                                 remove_cvref_t<alph_t>,
                                                 std::type_identity<alph_t>>>
struct char_is_valid_for_fn
{
private:
    SEQAN3_CPO_IMPL(3, (char_is_valid_for(v, s_alph_t{})                                       ))    // ADL
    SEQAN3_CPO_IMPL(2, (seqan3::adaptation::char_is_valid_for(v, s_alph_t{})                   ))    // customisation ns
    SEQAN3_CPO_IMPL(1, (deferred_type_t<remove_cvref_t<alph_t>, decltype(v)>::char_is_valid(v) ))    // member
    SEQAN3_CPO_IMPL(0, (to_char(assign_char_to(v, s_alph_t{})) == v                            ))    // fallback

public:
    //!\brief Operator definition.
    template <typename dummy = int> // need to make this a template to enforce deferred instantiation
    constexpr bool operator()(alphabet_char_t<alph_t> const a) const noexcept
    //!\cond
        requires requires (alphabet_char_t<alph_t> const a) { { impl(priority_tag<3>{}, a, dummy{}) }; }
    //!\endcond
    {
        static_assert(noexcept(impl(priority_tag<3>{}, a)),
            "Only overloads that are marked noexcept are picked up by seqan3::char_is_valid_for().");
        static_assert(std::Same<bool, decltype(impl(priority_tag<3>{}, a))>,
            "The return type of your char_is_valid_for() implementation must be 'bool'.");

        return impl(priority_tag<3>{}, a);
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Returns whether a character is in the valid set of a seqan3::Alphabet (usually implies a bijective mapping
 *        to an alphabet value).
 * \tparam your_type The alphabet type being queried.
 * \param  chr       The character being checked; must be convertible to `seqan3::alphabet_char_t<your_type>`.
 * \param  alph      The target object; its type must model seqan3::Alphabet.
 * \returns `true` or `false`.
 * \ingroup alphabet
 *
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A free function `char_is_valid_for(char_type const chr, your_type const &)` in the namespace of your
 *      type (or as `friend`). The function must be marked `noexcept` (`constexpr` is not required,
 *      but recommended) and the return type be `bool`. The value of the second argument to the function shall be
 *      ignored, it is only used to select the function via
 *      [argument-dependent lookup](https://en.cppreference.com/w/cpp/language/adl).
 *   2. A free function `char_is_valid_for(char_type const chr, your_type const &)` in
 *      `namespace seqan3::adaptation`. Used only by the SeqAn library internally. The same restrictions apply as above.
 *   3. A `static` member function called `char_is_valid(char_type)` (not `char_is_valid_for`). It must
 *      be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be `bool`.
 *
 * An alphabet type *may* provide one of the above. If none is provided, this function will declare every character
 * `c` as valid for whom it holds that `seqan3::to_char(seqan3::assign_char_to(c, alph_t{})) == c`, i.e. converting
 * back and forth results in the same value.
 *
 * *Note* that if the alphabet type with cvref removed is not std::is_nothrow_default_constructible, this function
 * object will instead look for `char_is_valid_for(char_type const chr, std::type_identity<your_type> const &)`
 * with the same semantics. In that case the "fallback" above also does not work and you are required to provide
 * such an implementation.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/char_is_valid_for.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (TODO link to manual). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
template <typename alph_t>
//!\cond
    requires requires { { to_char(std::declval<alph_t>()) }; } // to_char() is required by some defs
//!\endcond
inline constexpr auto char_is_valid_for = detail::adl::only::char_is_valid_for_fn<alph_t>{};
//!\}
} // namespace seqan3

// ============================================================================
// assign_char_strictly_to()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::assign_char_strictly_to.
//!\ingroup alphabet
struct assign_char_strictly_to_fn
{
    //!\brief Operator overload for lvalues.
    template <typename alph_t>
    //!\cond
        requires requires (alph_t a, seqan3::alphabet_char_t<alph_t> r)
        {
            { seqan3::assign_char_to(r, a) } -> alph_t;
            { seqan3::char_is_valid_for<alph_t>(r) } -> bool;
        }
    //!\endcond
    decltype(auto) operator()(seqan3::alphabet_char_t<alph_t> const r, alph_t & a) const
    {
        if (!seqan3::char_is_valid_for<alph_t>(r))
            //TODO: instead of doing a narrowing cast here, add constructor overload for other chars?
            throw seqan3::invalid_char_assignment{seqan3::detail::get_display_name_v<alph_t>, static_cast<char>(r)};

        return seqan3::assign_char_to(r, a);
    }

    //!\brief Operator overload for rvalues.
    template <typename alph_t>
    //!\cond
        requires requires (alph_t a, seqan3::alphabet_char_t<alph_t> r)
        {
            { seqan3::assign_char_to(r, a) } -> alph_t;
            { seqan3::char_is_valid_for<alph_t>(r) } -> bool;
        }
    //!\endcond
    auto operator()(seqan3::alphabet_char_t<alph_t> const r, alph_t && a) const
    {
        return operator()(r, a); // call above function but return by value
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects
 * \{
 */

/*!\brief Assign a character to an alphabet object, throw if the character is not valid.
 * \tparam your_type Type of the target object.
 * \param chr  The character being assigned; must be of the seqan3::alphabet_char_t of the target object.
 * \param alph The target object; its type must model seqan3::Alphabet.
 * \returns Reference to `alph` if `alph` was given as lvalue, otherwise a copy.
 * \throws seqan3::invalid_char_assignment If `seqan3::char_is_valid_for<decltype(alph)>(chr) == false`.
 * \ingroup alphabet

 * \details
 *
 * This is a function object. Invoke it with the parameters specified above.
 *
 * Note that this is not a customisation point and it cannot be "overloaded".
 * It simply invokes seqan3::char_is_valid_for and seqan3::assign_char_to.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/assign_char_strictly_to.cpp
 *
 */
inline constexpr auto assign_char_strictly_to = detail::adl::only::assign_char_strictly_to_fn{};
//!\}
} // namespace seqan3

// ============================================================================
// alphabet_size_v
// ============================================================================

namespace seqan3::detail::adl::only
{
/*!\brief Functor definition used indirectly by for seqan3::detail::alphabet_size_v.
 * \tparam alph_t   The type being queried.
 * \tparam s_alph_t `alph_t` with cvref removed and possibly wrapped in std::type_identity; never user-provide this!
 * \ingroup alphabet
 */
template <typename alph_t,
          typename s_alph_t = std::conditional_t<std::is_nothrow_default_constructible_v<remove_cvref_t<alph_t>> &&
                                                 seqan3::is_constexpr_default_constructible_v<remove_cvref_t<alph_t>>,
                                                 remove_cvref_t<alph_t>,
                                                 std::type_identity<alph_t>>>
struct alphabet_size_fn
{
private:
    SEQAN3_CPO_IMPL(2, (alphabet_size(v)                           ))    // ADL
    SEQAN3_CPO_IMPL(1, (seqan3::adaptation::alphabet_size(v)       ))    // customisation namespace
    SEQAN3_CPO_IMPL(0, (deferred_type_t<remove_cvref_t<alph_t>, decltype(v)>::alphabet_size )) // member

public:
    //!\brief Operator definition.
    template <typename dummy = int> // need to make this a template to enforce deferred instantiation
    constexpr auto operator()() const noexcept
    //!\cond
        requires requires { { impl(priority_tag<2>{}, s_alph_t{}, dummy{}) }; }
    //!\endcond
    {
        static_assert(noexcept(impl(priority_tag<2>{}, s_alph_t{})),
            "Only overloads that are marked noexcept are picked up by seqan3::alphabet_size.");
        static_assert(std::Constructible<size_t, decltype(impl(priority_tag<2>{}, s_alph_t{}))>,
            "The return type of your alphabet_size implementation must be convertible to size_t.");
        static_assert(SEQAN3_IS_CONSTEXPR(impl(priority_tag<2>{}, s_alph_t{})),
            "Only overloads that are marked constexpr are picked up by seqan3::alphabet_size.");

        return impl(priority_tag<2>{}, s_alph_t{});
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\brief A type trait that holds the size of a (semi-)alphabet.
 * \tparam your_type The (semi-)alphabet type being queried.
 * \ingroup alphabet
 *
 * \details
 *
 * This type trait is implemented as a global variable template.
 *
 * It is only defined for types that provide one of the following (checked in this order):
 *
 *   1. A free function `alphabet_size(your_type const &)` in the namespace of your type (or as `friend`) that
 *      returns the size as an integral value. The function must be marked `constexpr` and `noexcept` and the return
 *      type needs to be implicitly convertible to `size_t`. The value of the argument to the function shall be ignored,
 *      it is only used to select the function via
 *      [argument-dependent lookup](https://en.cppreference.com/w/cpp/language/adl).
 *   2. A free function `alphabet_size(your_type const &)` in `namespace seqan3::adaptation` that returns
 *      the size as an integral value. Used only by the SeqAn library internally. The same restrictions apply as above.
 *   3. A `static constexpr` data member called `alphabet_size` that is the size. It must
 *      be implicitly convertible to `size_t`.
 *
 * Every (semi-)alphabet type must provide one of the above.
 *
 * *Note* that if the (semi-)alphabet type with cvref removed is not std::is_nothrow_default_constructible or not
 * seqan3::is_constexpr_default_constructible, this object will instead look for
 * `alphabet_size(std::type_identity<your_type> const &)` with the same semantics (in cases 1. and 2.).
 *
 * ### Example
 *
 * \include test/snippet/alphabet/alphabet_size.cpp
 *
 * For an example of a full alphabet definition with free function implementations (solution 1. above),
 * see seqan3::Alphabet.
 *
 * ### Customisation point
 *
 * This is a customisation point (TODO link to manual). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
template <typename alph_t>
//!\cond
    requires requires { { detail::adl::only::alphabet_size_fn<alph_t>{}() }; }
//!\endcond
inline constexpr auto alphabet_size_v = detail::adl::only::alphabet_size_fn<alph_t>{}();

// ============================================================================
// Semialphabet
// ============================================================================

/*!\interface seqan3::Semialphabet <>
 * \brief The basis for seqan3::Alphabet, but requires only rank interface (not char).
 * \extends std::StrictTotallyOrdered
 * \extends std::CopyConstructible
 * \ingroup alphabet
 *
 * This concept represents the "rank part" of what is considered "an alphabet" in SeqAn. It requires no
 * `char` representation and corresponding interfaces. It is mostly used internally.
 *
 * ### Requirements
 *
 *   1. `t` shall model std::StrictTotallyOrdered ("has all comparison operators")
 *   2. objects of type `t` shall be efficiently copyable:
 *     * `t` shall model std::CopyConstructible and be std::is_nothrow_copy_constructible
 *     * move construction shall not be more efficient than copy construction; this implies no dynamic memory
 *       (de-)allocation [this is a semantic requirement that cannot be checked]
 *   3. seqan3::alphabet_size_v needs to be defined for `t`
 *   4. seqan3::to_rank needs to be defined for objects of type `t`
 *
 * See the documentation pages for the respective requirements.
 * The implications of 2. are that you can always take function arguments of types that model seqan3::Semialphabet
 * by value.
 *
 * It is highly recommended that non-reference types that model this concept, also model:
 *
 *   * std::Regular
 *   * std::is_trivially_copyable
 *   * seqan3::StandardLayout
 *
 * All alphabets available in SeqAn (with very few exceptions) do so.
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
SEQAN3_CONCEPT Semialphabet =
    std::StrictTotallyOrdered<t> &&
    std::CopyConstructible<t> &&
    std::is_nothrow_copy_constructible_v<t> &&
    requires (t v)
{
    requires seqan3::alphabet_size_v<t> >= 0;

    { seqan3::to_rank(v) };
};
//!\endcond

// ============================================================================
// WritableSemialphabet
// ============================================================================

/*!\interface seqan3::WritableSemialphabet <>
 * \brief A refinement of seqan3::Semialphabet that adds assignability.
 * \extends seqan3::Semialphabet
 * \ingroup alphabet
 *
 * This concept refines seqan3::Semialphabet and adds the requirement to be able to change the value by
 * assigning a value of the rank representation.
 *
 * For a detailed overview of how the different alphabet concepts are related, see
 * \ref alphabet module.
 *
 * ### Requirements
 *
 *   1. `t` shall model seqan3::Semialphabet
 *   2. seqan3::assign_rank_to needs to be defined for objects of type `t`
 *
 * See the documentation pages for the respective requirements.
 *
 * ### Related types
 *
 * If a given type `t` models this concept, the following types typically do so, as well:
 *
 *   * `t &`
 *
 * `const`-qualified types on the other hand are not assignable.
 *
 * ### Serialisation
 *
 * Types that model the concept (and all refinements) can be serialised via SeqAn3
 * serialisation support.
 * The rank value is (de-)serialised, types need not provide any overloads themselves.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT WritableSemialphabet = Semialphabet<t> && requires (t v)
{
    { seqan3::assign_rank_to(0, v) };
};
//!\endcond

// ============================================================================
// Alphabet
// ============================================================================

/*!\interface seqan3::Alphabet <>
 * \brief The generic alphabet concept that covers most data types used in ranges.
 * \extends seqan3::Semialphabet
 * \ingroup alphabet
 *
 * This is the core alphabet concept that many other alphabet concepts refine.
 *
 * For a detailed overview of how the different alphabet concepts are related, see
 * \ref alphabet module.
 *
 * ### Requirements
 *
 *   1. `t` shall model seqan3::Semialphabet ("has all rank representation")
 *   4. seqan3::to_char needs to be defined for objects of type `t`
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
 * ### Writing your own alphabet
 *
 * This is an example of a minimal custom alphabet that provides implementations for all necessary customisation
 * points:
 *
 * \snippet test/unit/alphabet/custom_alphabet_test.cpp my_alph
 *
 * This is an example of a custom alphabet that is not default-constructible and that has a non-default overload for
 * seqan3::char_is_valid_for:
 *
 * \snippet test/unit/alphabet/custom_alphabet2_test.cpp my_alph
 *
 * (Note that you should really make your alphabet types no-throw-default-constructible if you can!)
 *
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT Alphabet = Semialphabet<t> && requires (t v)
{
    { seqan3::to_char(v) };
};
//!\endcond

// ============================================================================
// WritableAlphabet
// ============================================================================

/*!\interface seqan3::WritableAlphabet <>
 * \brief Refines seqan3::Alphabet and adds assignability.
 * \extends seqan3::Alphabet
 * \extends seqan3::WritableSemialphabet
 * \ingroup alphabet
 *
 * This concept refines seqan3::Alphabet and seqan3::WritableSemialphabet and adds the requirement to be able to change
 * the value by assigning a value of the character representation.
 *
 * For a detailed overview of how the different alphabet concepts are related, see
 * \ref alphabet module.
 *
 * ### Requirements
 *
 *   1. `t` shall model seqan3::Alphabet
 *   2. `t` shall model seqan3::WritableSemialphabet
 *   3. seqan3::assign_char_to needs to be defined for objects of type `t`
 *   4. seqan3::char_is_valid_for needs to be defined for type `t` and an argument of the character representation
 *
 * See the documentation pages for the respective requirements.
 *
 * ### Related types
 *
 * If a given type `t` models this concept, the following types typically do so, as well:
 *
 *   * `t &`
 *
 * `const`-qualified types on the other hand are not assignable.
 *
 * ### Serialisation
 *
 * Types that model the concept (and all refinements) can be serialised via SeqAn3
 * serialisation support.
 * The rank value is (de-)serialised, types need not provide any overloads themselves.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT WritableAlphabet = Alphabet<t> && WritableSemialphabet<t> && requires (t v)
{
    { seqan3::assign_char_to(0, v) };

    { seqan3::char_is_valid_for<t>(std::declval<decltype(seqan3::to_char(v))>()) };
};
//!\endcond

// ============================================================================
//  serialisation
// ============================================================================

/*!\cond DEV
 * \name Generic serialisation functions for all seqan3::Semialphabet
 * \brief All types that satisfy seqan3::Semialphabet can be serialised via Cereal.
 *
 * \{
 */
/*!
 * \brief Save an alphabet letter to stream.
 * \tparam archive_t Must satisfy seqan3::CerealOutputArchive.
 * \tparam alphabet_t Type of l; must satisfy seqan3::Semialphabet.
 * \param l The alphabet letter.
 * \relates seqan3::Semialphabet
 *
 * \details
 *
 * Delegates to seqan3::Semialphabet::to_rank.
 *
 * \attention These functions are never called directly, see the \ref alphabet module on how to use serialisation.
 */
template <CerealOutputArchive archive_t, Semialphabet alphabet_t>
alphabet_rank_t<alphabet_t> CEREAL_SAVE_MINIMAL_FUNCTION_NAME(archive_t const &, alphabet_t const & l)
{
    return to_rank(l);
}

/*!\brief Restore an alphabet letter from a saved rank.
 * \tparam archive_t Must satisfy seqan3::CerealInputArchive.
 * \tparam wrapped_alphabet_t A seqan3::Semialphabet after Cereal mangles it up.
 * \param l The alphabet letter (cereal wrapped).
 * \param r The assigned value.
 * \relates seqan3::Semialphabet
 *
 * \details
 *
 * Delegates to seqan3::Semialphabet::assign_rank.
 *
 * \attention These functions are never called directly, see the \ref alphabet module on how to use serialisation.
 */
template <CerealInputArchive archive_t, typename wrapped_alphabet_t>
void CEREAL_LOAD_MINIMAL_FUNCTION_NAME(archive_t const &,
                                       wrapped_alphabet_t && l,
                                       alphabet_rank_t<detail::strip_cereal_wrapper_t<wrapped_alphabet_t>> const & r)
    requires Semialphabet<detail::strip_cereal_wrapper_t<wrapped_alphabet_t>>
{
    assign_rank_to(r, static_cast<detail::strip_cereal_wrapper_t<wrapped_alphabet_t> &>(l));
}
/*!\}
 * \endcond
 */

} // namespace seqan3

namespace seqan3::detail
{
// ============================================================================
// ConstexprSemialphabet
// ============================================================================

/*!\interface seqan3::detail::ConstexprSemialphabet <>
 * \brief A seqan3::Semialphabet that has constexpr accessors.
 * \extends seqan3::Semialphabet
 * \ingroup alphabet
 *
 * The same as seqan3::Semialphabet, except that all required functions are also required to be callable
 * in a `constexpr`-context.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT ConstexprSemialphabet = Semialphabet<t> && requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(to_rank(std::remove_reference_t<t>{}));
};
//!\endcond

// ============================================================================
// WritableConstexprSemialphabet
// ============================================================================

/*!\interface seqan3::detail::WritableConstexprSemialphabet <>
 * \brief A seqan3::WritableSemialphabet that has a constexpr assignment.
 * \extends seqan3::detail::ConstexprSemialphabet
 * \extends seqan3::WritableSemialphabet
 * \ingroup alphabet
 *
 * Refines seqan3::detail::ConstexprSemialphabet and seqan3::WritableSemialphabet and requires that the call to
 * seqan3::assign_rank_to can happen in a `constexpr`-context.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT WritableConstexprSemialphabet = ConstexprSemialphabet<t> && WritableSemialphabet<t> && requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(seqan3::assign_rank_to(0, std::remove_reference_t<t>{}));
};
//!\endcond

// ============================================================================
// ConstexprAlphabet
// ============================================================================

/*!\interface seqan3::detail::ConstexprAlphabet <>
 * \brief A seqan3::Alphabet that has constexpr accessors.
 * \extends seqan3::detail::ConstexprSemialphabet
 * \extends seqan3::Alphabet
 * \ingroup alphabet
 *
 * Refines seqan3::detail::ConstexprSemialphabet and seqan3::Alphabet and requires that the call to
 * seqan3::to_char can happen in a `constexpr`-context.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT ConstexprAlphabet = ConstexprSemialphabet<t> && Alphabet<t> && requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(to_char(std::remove_reference_t<t>{}));
};
//!\endcond

// ============================================================================
// WritableConstexprAlphabet
// ============================================================================

/*!\interface seqan3::detail::WritableConstexprAlphabet <>
 * \brief A seqan3::WritableAlphabet that has constexpr accessors.
 * \extends seqan3::detail::ConstexprAlphabet
 * \extends seqan3::detail::WritableConstexprSemialphabet
 * \extends seqan3::WritableAlphabet
 * \ingroup alphabet
 *
 * Refines seqan3::detail::ConstexprAlphabet, seqan3::detail::WritableConstexprSemialphabet and
 * seqan3::WritableAlphabet and requires that the calls to seqan3::assign_char_to and seqan3::char_is_valid_for
 * can happen in a `constexpr`-context.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT WritableConstexprAlphabet =
    ConstexprAlphabet<t> && WritableConstexprSemialphabet<t> && WritableAlphabet<t> && requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(seqan3::assign_char_to(alphabet_char_t<t>{}, std::remove_reference_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(seqan3::char_is_valid_for<t>(alphabet_char_t<t>{}));
};
//!\endcond

} // namespace seqan3::detail

#include <seqan3/alphabet/detail/hash.hpp>
