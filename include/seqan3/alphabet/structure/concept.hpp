// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::RnaStructureAlphabet.
 */

#pragma once

#include <optional>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/customisation_point.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

// ============================================================================
// forwards
// ============================================================================

//!\cond
namespace seqan3::custom
{

void is_pair_open();
void is_pair_close();
void is_unpaired();
void max_pseudoknot_depth();
void pseudoknot_id();

} // namespace seqan3::custom
//!\endcond

// ============================================================================
// is_pair_open()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::is_pair_open.
struct is_pair_open_fn
{
private:
    SEQAN3_CPO_IMPL(2, is_pair_open(v)                     ) // ADL
    SEQAN3_CPO_IMPL(1, seqan3::custom::is_pair_open(v)     ) // customisation namespace
    SEQAN3_CPO_IMPL(0, v.is_pair_open()                    ) // member

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
        static_assert(std::Same<bool, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your is_pair_open() implementation must be 'bool'.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl::only

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
 *   1. A free function `is_pair_open(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be `bool`.
 *   2. A free function `is_pair_open(your_type const a)` in `namespace seqan3::custom`.
 *      The same restrictions apply as above.
 *   3. A member function called `is_pair_open()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `bool`.
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
inline constexpr auto is_pair_open = detail::adl::only::is_pair_open_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// is_pair_close()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::is_pair_close.
struct is_pair_close_fn
{
private:
    SEQAN3_CPO_IMPL(2, is_pair_close(v)                     ) // ADL
    SEQAN3_CPO_IMPL(1, seqan3::custom::is_pair_close(v)     ) // customisation namespace
    SEQAN3_CPO_IMPL(0, v.is_pair_close()                    ) // member

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
        static_assert(std::Same<bool, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your is_pair_close() implementation must be 'bool'.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl::only

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
 *   1. A free function `is_pair_close(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be `bool`.
 *   2. A free function `is_pair_close(your_type const a)` in `namespace seqan3::custom`.
 *      The same restrictions apply as above.
 *   3. A member function called `is_pair_close()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `bool`.
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
inline constexpr auto is_pair_close = detail::adl::only::is_pair_close_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// is_unpaired()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::is_unpaired.
struct is_unpaired_fn
{
private:
    SEQAN3_CPO_IMPL(2, is_unpaired(v)                     ) // ADL
    SEQAN3_CPO_IMPL(1, seqan3::custom::is_unpaired(v)     ) // customisation namespace
    SEQAN3_CPO_IMPL(0, v.is_unpaired()                    ) // member

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
        static_assert(std::Same<bool, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your is_unpaired() implementation must be 'bool'.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl::only

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
 *   1. A free function `is_unpaired(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be `bool`.
 *   2. A free function `is_unpaired(your_type const a)` in `namespace seqan3::custom`.
 *      The same restrictions apply as above.
 *   3. A member function called `is_unpaired()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `bool`.
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
inline constexpr auto is_unpaired = detail::adl::only::is_unpaired_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// max_pseudoknot_depth
// ============================================================================

namespace seqan3::detail::adl::only
{

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
private:
    SEQAN3_CPO_IMPL(2, (max_pseudoknot_depth(v)                                                    )) // ADL
    SEQAN3_CPO_IMPL(1, (seqan3::custom::max_pseudoknot_depth(v)                                    )) // custom nsp
    SEQAN3_CPO_IMPL(0, (deferred_type_t<remove_cvref_t<alph_t>, decltype(v)>::max_pseudoknot_depth )) // member

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
        static_assert(std::Constructible<size_t, decltype(impl(priority_tag<2>{}, s_alph_t{}))>,
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

} // namespace seqan3::detail::adl::only

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
 *   1. A free function `max_pseudoknot_depth(your_type const)` in the namespace of your type (or as `friend`).
 *      The function must be marked `constexpr` and `noexcept` and the return type must be convertible to `size_t`.
 *      The value of the argument to the function shall be ignored, it is only used to select the function via
 *      [argument-dependent lookup](https://en.cppreference.com/w/cpp/language/adl).
 *   2. A free function `max_pseudoknot_depth(your_type const)` in `namespace seqan3::custom`.
 *      The same restrictions apply as above.
 *   3. A `static constexpr` data member of a type implicitly convertible to `size_t` called `max_pseudoknot_depth`.
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
    requires requires { { detail::adl::only::max_pseudoknot_depth_fn<alph_t>{} }; } &&
             requires { { detail::adl::only::max_pseudoknot_depth_obj<alph_t>() }; }
//!\endcond
inline constexpr auto max_pseudoknot_depth = detail::adl::only::max_pseudoknot_depth_obj<alph_t>();

} // namespace seqan3

// ============================================================================
// pseudoknot_id()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::pseudoknot_id.
struct pseudoknot_id_fn
{
private:
    SEQAN3_CPO_IMPL(2, pseudoknot_id(v)                     ) // ADL
    SEQAN3_CPO_IMPL(1, seqan3::custom::pseudoknot_id(v)     ) // customisation namespace
    SEQAN3_CPO_IMPL(0, v.pseudoknot_id()                    ) // member

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
        static_assert(std::Constructible<std::optional<size_t>, decltype(impl(priority_tag<2>{}, chr))>,
            "The return type of your pseudoknot_id() implementation must be convertible to std::optional<size_t>.");

        return impl(priority_tag<2>{}, chr);
    }
};

} // namespace seqan3::detail::adl::only

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
 *   1. A free function `pseudoknot_id(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type must be convertible to `size_t`.
 *   2. A free function `pseudoknot_id(your_type const a)` in `namespace seqan3::custom`.
 *      The same restrictions apply as above.
 *   3. A member function called `pseudoknot_id()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type
 *      must be convertible to `size_t`.
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
inline constexpr auto pseudoknot_id = detail::adl::only::pseudoknot_id_fn{};
//!\}

} // namespace seqan3

// ============================================================================
// RnaStructureAlphabet concept
// ============================================================================

namespace seqan3
{
/*!\interface seqan3::RnaStructureAlphabet <>
 * \brief A concept that indicates whether an alphabet represents RNA structure.
 * \extends seqan3::Alphabet
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
 *   1. `t` shall model seqan3::Alphabet
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
SEQAN3_CONCEPT RnaStructureAlphabet = seqan3::Alphabet<t> && requires(t val)
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
