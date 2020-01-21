// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::rna_structure_alphabet.
 */

#pragma once

#include <optional>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/customisation_point.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

// ============================================================================
// is_pair_open()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename ...args_t>
void is_pair_open(args_t ...) = delete;

//!\brief Functor definition for seqan3::is_pair_open.
struct is_pair_open_fn
{
public:
    SEQAN3_CPO_IMPL(2, seqan3::custom::alphabet<decltype(v)>::is_pair_open(v)) // explicit customisation
    SEQAN3_CPO_IMPL(1, is_pair_open(v)                                       ) // ADL
    SEQAN3_CPO_IMPL(0, v.is_pair_open()                                      ) // member

public:
    //!\brief Operator definition.
    template <typename rna_structure_t>
    //!\cond
        requires requires (rna_structure_t const chr) { { impl(priority_tag<2>{}, chr) }; }
    //!\endcond
    constexpr auto operator()(rna_structure_t const chr) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, chr)),
            "Only overloads that are marked noexcept are picked up by seqan3::is_pair_open().");
        static_assert(std::same_as<bool, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your is_pair_open() implementation must be 'bool'.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\name Function objects (Structure)
 * \{
 */

/*!\brief Check whether the given character represents a rightward interaction in an RNA structure.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns True if the character represents a rightward interaction, False otherwise.
 * \ingroup structure
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A static member function `is_pair_open(your_type const a)` of the class `seqan3::custom::alphabet<your_type>`.
 *   2. A free function `is_pair_open(your_type const a)` in the namespace of your type (or as `friend`).
 *   3. A member function called `is_pair_open()`.
 *
 * Functions are only considered for one of the above cases if they are marked `noexcept` (`constexpr` is not required,
 * but recommended) and if the returned type is `bool.
 *
 * Every RNA structure alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/structure/wuss_is_pair_open.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto is_pair_open = detail::adl_only::is_pair_open_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// is_pair_close()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename ...args_t>
void is_pair_close(args_t ...) = delete;

//!\brief Functor definition for seqan3::is_pair_close.
struct is_pair_close_fn
{
public:
    SEQAN3_CPO_IMPL(2, seqan3::custom::alphabet<decltype(v)>::is_pair_close(v)) // explicit customisation
    SEQAN3_CPO_IMPL(1, is_pair_close(v)                                       ) // ADL
    SEQAN3_CPO_IMPL(0, v.is_pair_close()                                      ) // member

public:
    //!\brief Operator definition.
    template <typename rna_structure_t>
    //!\cond
        requires requires (rna_structure_t const chr) { { impl(priority_tag<2>{}, chr) }; }
    //!\endcond
    constexpr auto operator()(rna_structure_t const chr) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, chr)),
            "Only overloads that are marked noexcept are picked up by seqan3::is_pair_close().");
        static_assert(std::same_as<bool, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your is_pair_close() implementation must be 'bool'.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\name Function objects (Structure)
 * \{
 */

/*!\brief Check whether the given character represents a leftward interaction in an RNA structure.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns True if the character represents a leftward interaction, False otherwise.
 * \ingroup structure
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A static member function `is_pair_close(your_type const a)` of the class `seqan3::custom::alphabet<your_type>`.
 *   2. A free function `is_pair_close(your_type const a)` in the namespace of your type (or as `friend`).
 *   3. A member function called `is_pair_close()`.
 *
 * Functions are only considered for one of the above cases if they are marked `noexcept` (`constexpr` is not required,
 * but recommended) and if the returned type is `bool.
 *
 * Every RNA structure alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/structure/wuss_is_pair_close.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto is_pair_close = detail::adl_only::is_pair_close_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// is_unpaired()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename ...args_t>
void is_unpaired(args_t ...) = delete;

//!\brief Functor definition for seqan3::is_unpaired.
struct is_unpaired_fn
{
public:
    SEQAN3_CPO_IMPL(2, seqan3::custom::alphabet<decltype(v)>::is_unpaired(v)) // explicit customisation
    SEQAN3_CPO_IMPL(1, is_unpaired(v)                                       ) // ADL
    SEQAN3_CPO_IMPL(0, v.is_unpaired()                                      ) // member

public:
    //!\brief Operator definition.
    template <typename rna_structure_t>
    //!\cond
        requires requires (rna_structure_t const chr) { { impl(priority_tag<2>{}, chr) }; }
    //!\endcond
    constexpr auto operator()(rna_structure_t const chr) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, chr)),
            "Only overloads that are marked noexcept are picked up by seqan3::is_unpaired().");
        static_assert(std::same_as<bool, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your is_unpaired() implementation must be 'bool'.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\name Function objects (Structure)
 * \{
 */

/*!\brief Check whether the given character represents an unpaired nucleotide in an RNA structure.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns True if the character represents an unpaired site, False otherwise.
 * \ingroup structure
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A static member function `is_unpaired(your_type const a)` of the class `seqan3::custom::alphabet<your_type>`.
 *   2. A free function `is_unpaired(your_type const a)` in the namespace of your type (or as `friend`).
 *   3. A member function called `is_unpaired()`.
 *
 * Functions are only considered for one of the above cases if they are marked `noexcept` (`constexpr` is not required,
 * but recommended) and if the returned type is `bool.
 *
 * Every RNA structure alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/structure/wuss_is_unpaired.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto is_unpaired = detail::adl_only::is_unpaired_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// max_pseudoknot_depth
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename ...args_t>
void max_pseudoknot_depth(args_t ...) = delete;

/*!\brief Functor definition for seqan3::max_pseudoknot_depth.
 * \tparam alph_t   The type being queried.
 * \tparam s_alph_t `alph_t` with cvref removed and possibly wrapped in std::type_identity; never user-provide this!
 * \ingroup structure
 */
template <typename alph_t,
          typename s_alph_t = std::conditional_t<std::is_nothrow_default_constructible_v<remove_cvref_t<alph_t>> &&
                                                 seqan3::is_constexpr_default_constructible_v<remove_cvref_t<alph_t>>,
                                                 remove_cvref_t<alph_t>,
                                                 std::type_identity<alph_t>>>
struct max_pseudoknot_depth_fn
{
public:
    SEQAN3_CPO_IMPL(2, (deferred_type_t<seqan3::custom::alphabet<alph_t>, decltype(v)>::max_pseudoknot_depth)) // custom
    SEQAN3_CPO_IMPL(1, (max_pseudoknot_depth(v)                                                             )) // ADL
    SEQAN3_CPO_IMPL(0, (deferred_type_t<remove_cvref_t<alph_t>, decltype(v)>::max_pseudoknot_depth          )) // member

public:
    //!\brief Operator definition.
    template <typename dummy = int>
    //!\cond
        requires requires { { impl(priority_tag<2>{}, s_alph_t{}, dummy{}) }; }
    //!\endcond
    constexpr auto operator()() const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, s_alph_t{})),
            "Only overloads that are marked noexcept are picked up by seqan3::max_pseudoknot_depth.");
        static_assert(std::constructible_from<size_t, decltype(impl(priority_tag<2>{}, s_alph_t{}))>,
            "The return type of your max_pseudoknot_depth implementation must be convertible to size_t.");
        static_assert(SEQAN3_IS_CONSTEXPR(impl(priority_tag<2>{}, s_alph_t{})),
            "Only overloads that are marked constexpr are picked up by seqan3::max_pseudoknot_depth.");

        return impl(priority_tag<2>{}, s_alph_t{});
    }
};

//!\cond
// required to prevent https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89953
template <typename alph_t>
    requires requires { { max_pseudoknot_depth_fn<alph_t>{} }; }
inline constexpr auto max_pseudoknot_depth_obj = max_pseudoknot_depth_fn<alph_t>{};
//!\endcond

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\name Function objects (Structure)
 * \{
 */

/*!\brief A type trait that holds the ability of the structure alphabet to represent pseudoknots,
 *        i.e. crossing interactions, up to a certain depth.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns The maximum supported nestedness, or 1 if the alphabet cannot support pseudoknots.
 * \ingroup structure
 * \details
 *
 * The value is the maximum allowed depth of pseudoknots.
 * A value of 1 denotes no pseudoknots `((....))`,
 * while higher values denote the maximum allowed complexity of
 * crossing interactions, e.g. depth 2 `(({....))}` or depth 3 `({[....)}]`.
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A static member variable `max_pseudoknot_depth` of the class `seqan3::custom::alphabet<your_type>`.
 *   2. A free function `max_pseudoknot_depth(your_type const)` in the namespace of your type (or as `friend`).
 *   3. A static member variable `max_pseudoknot_depth` of the class `your_type`.
 *
 * Functions/variables are only considered for one of the above cases if they are marked `noexcept` **and** `constexpr` and
 * if the returned type is convertible to `size_t`. For 2. the value of the argument to the function shall be
 * ignored, the argument is only used to select the function via
 * [argument-dependent lookup](https://en.cppreference.com/w/cpp/language/adl).
 *
 * Every RNA structure alphabet type must provide one of the above.
 *
 * ### Example
 *
 * These are the expressions to retrieve the value:
 * \include test/snippet/alphabet/structure/wuss_max_pseudoknot_depth.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
template <typename alph_t>
//!\cond
    requires requires { { detail::adl_only::max_pseudoknot_depth_fn<alph_t>{} }; } &&
             requires { { detail::adl_only::max_pseudoknot_depth_obj<alph_t>() }; }
//!\endcond
inline constexpr auto max_pseudoknot_depth = detail::adl_only::max_pseudoknot_depth_obj<alph_t>();

} // namespace seqan3

// ============================================================================
// pseudoknot_id()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename ...args_t>
void pseudoknot_id(args_t ...) = delete;

//!\brief Functor definition for seqan3::pseudoknot_id.
struct pseudoknot_id_fn
{
public:
    SEQAN3_CPO_IMPL(2, seqan3::custom::alphabet<decltype(v)>::pseudoknot_id(v)) // explicit customisation
    SEQAN3_CPO_IMPL(1, pseudoknot_id(v)                                       ) // ADL
    SEQAN3_CPO_IMPL(0, v.pseudoknot_id()                                      ) // member

public:
    //!\brief Operator definition.
    template <typename rna_structure_t>
    //!\cond
        requires requires (rna_structure_t const chr) { { impl(priority_tag<2>{}, chr) }; }
    //!\endcond
    constexpr auto operator()(rna_structure_t const chr) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, chr)),
            "Only overloads that are marked noexcept are picked up by seqan3::pseudoknot_id().");
        static_assert(std::constructible_from<std::optional<size_t>, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your pseudoknot_id() implementation must be convertible to std::optional<size_t>.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\name Function objects (Structure)
 * \{
 */

/*!\brief Retrieve an id for the level of a pseudoknotted interaction (also known as 'page number').
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns An std::optional containing the pseudoknot identifier if `alph` represents an interaction.
 * The returned value is std::nullopt for unpaired sites. For non-nested interactions the identifier is always 0.
 * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
 * \ingroup structure
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A static member function `pseudoknot_id(your_type const a)` of the class `seqan3::custom::alphabet<your_type>`.
 *   2. A free function `pseudoknot_id(your_type const a)` in the namespace of your type (or as `friend`).
 *   3. A member function of `your_type` called `pseudoknot_id()`.
 *
 * Functions are only considered for one of the above cases if they are marked `noexcept` (`constexpr` is not required,
 * but recommended) and the return type is convertible to `size_t`.
 *
 * Every RNA structure alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/structure/wuss_pseudoknot_id.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto pseudoknot_id = detail::adl_only::pseudoknot_id_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// rna_structure_alphabet concept
// ============================================================================

namespace seqan3
{
/*!\interface seqan3::rna_structure_alphabet <>
 * \brief A concept that indicates whether an alphabet represents RNA structure.
 * \extends seqan3::alphabet
 * \ingroup structure
 *
 * \details
 *
 * RNA structure alphabets are required to represent interactions among RNA nucleotides.
 * Therefore, each structure letter can be categorised as unpaired, opening an interaction, or closing an interaction.
 * Additionally, the ability of representing pseudoknots is a property of RNA structure types.
 *
 * ### Requirements
 *
 *   1. `t` shall model seqan3::alphabet
 *   2. seqan3::is_pair_open needs to be defined for objects of type `t`
 *   3. seqan3::is_pair_close needs to be defined for objects of type `t`
 *   4. seqan3::is_unpaired needs to be defined for objects of type `t`
 *   5. seqan3::max_pseudoknot_depth needs to be defined for `t` and be greater than zero
 *   6. seqan3::pseudoknot_id needs to be defined for objects of type `t`
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
SEQAN3_CONCEPT rna_structure_alphabet = seqan3::alphabet<t> && requires(t val)
{
    { seqan3::is_pair_open(val) };
    { seqan3::is_pair_close(val) };
    { seqan3::is_unpaired(val) };
    { seqan3::pseudoknot_id(val) };

    // this is delegated to a static class variable, which must not be 0
    requires seqan3::max_pseudoknot_depth<t> > 0;
};
//!\endcond

} // namespace seqan3
