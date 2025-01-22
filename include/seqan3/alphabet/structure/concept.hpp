// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::rna_structure_alphabet.
 */

#pragma once

#include <concepts>
#include <optional>
#include <type_traits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/customisation_point.hpp>

// ============================================================================
// is_pair_open()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void is_pair_open(args_t...) = delete;

/*!\brief seqan3::detail::customisation_point_object (CPO) definition for seqan3::is_pair_open.
 * \ingroup alphabet_structure
 */
struct is_pair_open_cpo : public detail::customisation_point_object<is_pair_open_cpo, 2>
{
    //!\brief CRTP base class seqan3::detail::customisation_point_object.
    using base_t = detail::customisation_point_object<is_pair_open_cpo, 2>;
    //!\brief Only this class is allowed to import the constructors from #base_t. (CRTP safety idiom)
    using base_t::base_t;

    /*!\brief CPO overload (check 1 out of 3): explicit customisation via `seqan3::custom::alphabet`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is a pair_open.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<2>, alphabet_t && alphabet)(
        /*return*/ seqan3::custom::alphabet<alphabet_t>::is_pair_open(std::forward<alphabet_t>(alphabet)) == true /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): argument dependent lookup (ADL), i.e. `is_pair_open(alphabet)`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is a pair_open.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<1>, alphabet_t && alphabet)(
        /*return*/ is_pair_open(std::forward<alphabet_t>(alphabet)) == true /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): member access, i.e. `alphabet.is_pair_open()`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is a pair_open.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<0>, alphabet_t && alphabet)(
        /*return*/ std::forward<alphabet_t>(alphabet).is_pair_open() == true /*;*/
    );
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{
/*!\brief Check whether the given character represents a rightward interaction in an RNA structure.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns True if the character represents a rightward interaction, False otherwise.
 * \ingroup alphabet_structure
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
 * but recommended) and if the returned type is `bool`.
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
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
inline constexpr auto is_pair_open = detail::adl_only::is_pair_open_cpo{};

} // namespace seqan3

// ============================================================================
// is_pair_close()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void is_pair_close(args_t...) = delete;

/*!\brief seqan3::detail::customisation_point_object (CPO) definition for seqan3::is_pair_close.
 * \ingroup alphabet_structure
 */
struct is_pair_close_cpo : public detail::customisation_point_object<is_pair_close_cpo, 2>
{
    //!\brief CRTP base class seqan3::detail::customisation_point_object.
    using base_t = detail::customisation_point_object<is_pair_close_cpo, 2>;
    //!\brief Only this class is allowed to import the constructors from #base_t. (CRTP safety idiom)
    using base_t::base_t;

    /*!\brief CPO overload (check 1 out of 3): explicit customisation via `seqan3::custom::alphabet`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is a pair_close.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<2>, alphabet_t && alphabet)(
        /*return*/ seqan3::custom::alphabet<alphabet_t>::is_pair_close(std::forward<alphabet_t>(alphabet)) == true /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): argument dependent lookup (ADL), i.e. `is_pair_close(alphabet)`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is a pair_close.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<1>, alphabet_t && alphabet)(
        /*return*/ is_pair_close(std::forward<alphabet_t>(alphabet)) == true /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): member access, i.e. `alphabet.is_pair_close()`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is a pair_close.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<0>, alphabet_t && alphabet)(
        /*return*/ std::forward<alphabet_t>(alphabet).is_pair_close() == true /*;*/
    );
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{
/*!\brief Check whether the given character represents a leftward interaction in an RNA structure.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns True if the character represents a leftward interaction, False otherwise.
 * \ingroup alphabet_structure
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
 * but recommended) and if the returned type is `bool`.
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
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
inline constexpr auto is_pair_close = detail::adl_only::is_pair_close_cpo{};

} // namespace seqan3

// ============================================================================
// is_unpaired()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void is_unpaired(args_t...) = delete;

/*!\brief seqan3::detail::customisation_point_object (CPO) definition for seqan3::is_unpaired.
 * \ingroup alphabet_structure
 */
struct is_unpaired_cpo : public detail::customisation_point_object<is_unpaired_cpo, 2>
{
    //!\brief CRTP base class seqan3::detail::customisation_point_object.
    using base_t = detail::customisation_point_object<is_unpaired_cpo, 2>;
    //!\brief Only this class is allowed to import the constructors from #base_t. (CRTP safety idiom)
    using base_t::base_t;

    /*!\brief CPO overload (check 1 out of 3): explicit customisation via `seqan3::custom::alphabet`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is unpaired.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<2>, alphabet_t && alphabet)(
        /*return*/ seqan3::custom::alphabet<alphabet_t>::is_unpaired(std::forward<alphabet_t>(alphabet)) == true /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): argument dependent lookup (ADL), i.e. `is_unpaired(alphabet)`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is unpaired.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<1>, alphabet_t && alphabet)(
        /*return*/ is_unpaired(std::forward<alphabet_t>(alphabet)) == true /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): member access, i.e. `alphabet.is_unpaired()`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet that is queried whether it is unpaired.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<0>, alphabet_t && alphabet)(
        /*return*/ std::forward<alphabet_t>(alphabet).is_unpaired() == true /*;*/
    );
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{
/*!\brief Check whether the given character represents an unpaired nucleotide in an RNA structure.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns True if the character represents an unpaired site, False otherwise.
 * \ingroup alphabet_structure
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
 * but recommended) and if the returned type is `bool`.
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
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
inline constexpr auto is_unpaired = detail::adl_only::is_unpaired_cpo{};

} // namespace seqan3

// ============================================================================
// max_pseudoknot_depth
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void max_pseudoknot_depth(args_t...) = delete;

/*!\brief seqan3::detail::customisation_point_object (CPO) definition for seqan3::max_pseudoknot_depth.
 * \tparam alphabet_t The alphabet type being queried.
 * \ingroup alphabet_structure
 */
template <typename alphabet_t>
struct max_pseudoknot_depth_cpo : public detail::customisation_point_object<max_pseudoknot_depth_cpo<alphabet_t>, 2>
{
    //!\brief CRTP base class seqan3::detail::customisation_point_object.
    using base_t = detail::customisation_point_object<max_pseudoknot_depth_cpo<alphabet_t>, 2>;
    //!\brief Only this class is allowed to import the constructors from #base_t. (CRTP safety idiom)
    using base_t::base_t;

    /*!\brief If `alphabet_type` isn't std::is_nothrow_default_constructible, max_pseudoknot_depth will be called with
     *        std::type_identity instead of a default constructed alphabet.
     */
    template <typename alphabet_type>
    using alphabet_or_type_identity =
        std::conditional_t<std::is_nothrow_default_constructible_v<std::remove_cvref_t<alphabet_type>>
                               && seqan3::is_constexpr_default_constructible_v<std::remove_cvref_t<alphabet_type>>,
                           std::remove_cvref_t<alphabet_type>,
                           std::type_identity<alphabet_type>>;

    /*!\brief CPO overload (check 1 out of 3): explicit customisation via `seqan3::custom::alphabet`
     * \tparam alphabet_type The type of the alphabet. (Needed to defer instantiation for incomplete types.)
     */
    template <typename alphabet_type = alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<2>)(
        /*return*/ seqan3::custom::alphabet<alphabet_type>::max_pseudoknot_depth /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): argument dependent lookup (ADL), i.e.
     *        `max_pseudoknot_depth(alphabet_type{})`
     * \tparam alphabet_type The type of the alphabet. (Needed to defer instantiation for incomplete types.)
     *
     * \details
     *
     * If the alphabet_type isn't std::is_nothrow_default_constructible,
     * `max_pseudoknot_depth(std::type_identity<alphabet_type>{})` will be called.
     */
    template <typename alphabet_type = alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<1>)(
        /*return*/ max_pseudoknot_depth(alphabet_or_type_identity<alphabet_type>{}) /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): static member access, i.e. `alphabet_type::max_pseudoknot_depth`
     * \tparam alphabet_type The type of the alphabet. (Needed to defer instantiation for incomplete types.)
     */
    template <typename alphabet_type = alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<0>)(
        /*return*/ std::remove_cvref_t<alphabet_type>::max_pseudoknot_depth /*;*/
    );
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{
/*!\brief A type trait that holds the ability of the structure alphabet to represent pseudoknots,
 *        i.e. crossing interactions, up to a certain depth.
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns The maximum supported nestedness, or 1 if the alphabet cannot support pseudoknots.
 * \ingroup alphabet_structure
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
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
template <typename alph_t>
    requires requires {
        { detail::adl_only::max_pseudoknot_depth_cpo<alph_t>{}() };
    }
inline constexpr auto max_pseudoknot_depth = detail::adl_only::max_pseudoknot_depth_cpo<alph_t>{}();

} // namespace seqan3

// ============================================================================
// pseudoknot_id()
// ============================================================================

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void pseudoknot_id(args_t...) = delete;

//!\brief seqan3::detail::customisation_point_object (CPO) definition for seqan3::pseudoknot_id_cpo.
//!\ingroup alphabet_structure
struct pseudoknot_id_cpo : public detail::customisation_point_object<pseudoknot_id_cpo, 2>
{
    //!\brief CRTP base class seqan3::detail::customisation_point_object.
    using base_t = detail::customisation_point_object<pseudoknot_id_cpo, 2>;
    //!\brief Only this class is allowed to import the constructors from #base_t. (CRTP safety idiom)
    using base_t::base_t;

    /*!\brief CPO overload (check 1 out of 3): explicit customisation via `seqan3::custom::alphabet`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet the pseudoknot_id is returned from.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<2>, alphabet_t && alphabet)(
        /*return*/ seqan3::custom::alphabet<alphabet_t>::pseudoknot_id(std::forward<alphabet_t>(alphabet)) /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): argument dependent lookup (ADL), i.e. `pseudoknot_id(alphabet)`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet the pseudoknot_id is returned from.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<1>, alphabet_t && alphabet)(
        /*return*/ pseudoknot_id(std::forward<alphabet_t>(alphabet)) /*;*/
    );

    /*!\brief CPO overload (check 1 out of 3): member access, i.e. `alphabet.pseudoknot_id()`
     * \tparam alphabet_t The type of the alphabet.
     * \param alphabet The alphabet the pseudoknot_id is returned from.
     */
    template <typename alphabet_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(priority_tag<0>, alphabet_t && alphabet)(
        /*return*/ std::forward<alphabet_t>(alphabet).pseudoknot_id() /*;*/
    );
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{
/*!\brief Retrieve an id for the level of a pseudoknotted interaction (also known as 'page number').
 * \tparam your_type Type of the argument.
 * \param  chr       The RNA structure character whose property is checked.
 * \returns A std::optional containing the pseudoknot identifier if `alph` represents an interaction.
 * The returned value is std::nullopt for unpaired sites. For non-nested interactions the identifier is always 0.
 * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
 * \ingroup alphabet_structure
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
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
inline constexpr auto pseudoknot_id = detail::adl_only::pseudoknot_id_cpo{};

} // namespace seqan3

// ============================================================================
// rna_structure_alphabet concept
// ============================================================================

namespace seqan3
{
/*!\interface seqan3::rna_structure_alphabet <>
 * \brief A concept that indicates whether an alphabet represents RNA structure.
 * \extends seqan3::alphabet
 * \ingroup alphabet_structure
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
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
//!\cond
template <typename t>
concept rna_structure_alphabet = seqan3::alphabet<t> && requires (t val) {
    { seqan3::is_pair_open(val) };
    { seqan3::is_pair_close(val) };
    { seqan3::is_unpaired(val) };
    { seqan3::pseudoknot_id(val) };

    // this is delegated to a static class variable, which must not be 0
    requires seqan3::max_pseudoknot_depth<t> > 0;
};
//!\endcond

} // namespace seqan3
