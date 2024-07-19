// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_tuple_base.
 */

#pragma once

#include <cassert>
#include <concepts>
#include <utility>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/composite/detail.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/detail/alphabet_proxy.hpp>
#include <seqan3/utility/detail/integer_traits.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/detail/type_list_algorithm.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/type_list/type_list.hpp>
#include <seqan3/utility/type_pack/traits.hpp>
#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>

namespace seqan3::detail
{

//!\brief Prevents wrong instantiations of seqan3::alphabet_tuple_base's equality comparison operators.
template <typename tuple_derived_t, typename rhs_t, typename... component_types>
inline constexpr bool tuple_general_guard =
    (!std::same_as<rhs_t, tuple_derived_t>) && (!std::same_as<rhs_t, alphabet_tuple_base<component_types...>>)
    && (!std::is_base_of_v<tuple_derived_t, rhs_t>) && (!(std::same_as<rhs_t, component_types> || ...))
    && (!list_traits::contains<tuple_derived_t, recursive_required_types_t<rhs_t>>);

//!\brief Prevents wrong instantiations of seqan3::alphabet_tuple_base's equality comparison operators.
template <typename lhs_t, typename tuple_derived_t, typename rhs_t, typename... component_types>
inline constexpr bool tuple_eq_guard =
    (instantiate_if_v<lazy<weakly_equality_comparable_with_trait, rhs_t, component_types>,
                      std::same_as<lhs_t, tuple_derived_t>
                          && tuple_general_guard<tuple_derived_t, rhs_t, component_types...>>
     || ...);

//!\brief Prevents wrong instantiations of seqan3::alphabet_tuple_base's ordered comparison operators.
template <typename lhs_t, typename tuple_derived_t, typename rhs_t, typename... component_types>
inline constexpr bool tuple_order_guard =
    (instantiate_if_v<lazy<weakly_ordered_with_trait, rhs_t, component_types>,
                      std::same_as<lhs_t, tuple_derived_t>
                          && tuple_general_guard<lhs_t, tuple_derived_t, rhs_t, component_types...>>
     || ...);

} // namespace seqan3::detail

namespace seqan3
{

// forwards
//!\cond
template <typename t>
decltype(auto) get();

template <size_t i>
decltype(auto) get();
//!\endcond

/*!\brief The CRTP base for a combined alphabet that contains multiple values of different alphabets at the same time.
 * \ingroup alphabet_composite
 * \implements seqan3::writable_semialphabet
 * \if DEV
 * \implements seqan3::detail::writable_constexpr_semialphabet
 * \tparam component_types Types of letters; must model seqan3::detail::writable_constexpr_semialphabet.
 * \else
 * \tparam component_types Types of letters; must model std::regular and seqan3::writable_semialphabet and all
 * required function calls need to be callable in `constexpr`-context.
 * \endif
 *
 *
 * This data structure is a CRTP base class for combined alphabets, where the different
 * alphabet letters exist independently as a components, similar to a tuple.
 *
 * Short description:
 *   * combines multiple alphabets as independent components, similar to a tuple;
 *   * models seqan3::tuple_like, i.e. provides a get interface to its component_list;
 *   * is itself a seqan3::writable_semialphabet, but most derived types implement the full seqan3::writable_alphabet;
 *   * its alphabet size is the product of the individual sizes;
 *   * constructible, assignable and comparable with each component type and also all types that
 *     these are constructible/assignable/comparable with;
 *   * explicitly convertible to each of its component types
 *
 * \attention
 * This is a "pure base class", you cannot instantiate it, you can only inherit from it.
 * Most likely you are interested in using one of it's descendants like seqan3::qualified or seqan3::masked.
 * \cond DEV
 * To make a derived class "complete", you should add at least the following:
 *   * .to_char() member
 *   * .assign_char() member
 *   * a type deduction guide
 * \endcond
 *
 * \sa qualified
 * \sa masked
 *
 * \stableapi{Since version 3.1.}
 */
template <typename derived_type, typename... component_types>
    requires (detail::writable_constexpr_semialphabet<component_types> && ...) && (std::regular<component_types> && ...)
class alphabet_tuple_base :
    public alphabet_base<derived_type,
                         (1 * ... * alphabet_size<component_types>),
                         void> // no char type, because this is only semi_alphabet
{
private:
    //!\brief The base type of this class.
    using base_t = alphabet_base<derived_type,
                                 (1 * ... * alphabet_size<component_types>),
                                 void>; // no char type, because this is only semi_alphabet

    //!\brief A seqan3::type_list The types of each component in the composite
    using component_list = seqan3::type_list<component_types...>;

    //!\brief Is set to `true` if the type is contained in the type list.
    template <typename type>
    static constexpr bool is_component = seqan3::list_traits::contains<type, component_list>;

    //!\brief Is set to `true` if the type is uniquely contained in the type list.
    template <typename type>
    static constexpr bool is_unique_component = (seqan3::list_traits::count<type, component_list> == 1);

    // forward declaration: see implementation below
    template <typename alphabet_type, size_t index>
    class component_proxy;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_tuple_base() noexcept : base_t{}
    {}                                                                                //!< Defaulted.
    constexpr alphabet_tuple_base(alphabet_tuple_base const &) = default;             //!< Defaulted.
    constexpr alphabet_tuple_base(alphabet_tuple_base &&) = default;                  //!< Defaulted.
    constexpr alphabet_tuple_base & operator=(alphabet_tuple_base const &) = default; //!< Defaulted.
    constexpr alphabet_tuple_base & operator=(alphabet_tuple_base &&) = default;      //!< Defaulted.
    ~alphabet_tuple_base() = default;                                                 //!< Defaulted.

    using base_t::base_t;
    //!\}

    //!\brief Befriend the derived type so that it can instantiate.
    //!\sa https://isocpp.org/blog/2017/04/quick-q-prevent-user-from-derive-from-incorrect-crtp-base
    friend derived_type;

    // Import from base:
    using typename base_t::rank_type;

public:
    // Import from base:
    using base_t::alphabet_size;
    using base_t::assign_rank;
    using base_t::to_rank;

    //!\brief Export this type's components in a visible manner.
    //!\private
    using seqan3_required_types = component_list;
    //!\brief Export this type's components and possibly the components' components in a visible manner.
    //!\private
    using seqan3_recursive_required_types = list_traits::concat<
        component_list,
        detail::transformation_trait_or_t<detail::recursive_required_types<component_types>, seqan3::type_list<>>...>;
    //!\brief Make specialisations of this template identifiable in metapgrogramming contexts.
    //!\private
    static constexpr bool seqan3_alphabet_tuple_like = true;

    /*!\name Constructors, destructor and assignment
     * \{
     * \attention Please do not directly use the CRTP base class. The functions
     *            are only public for the usage in their derived classes (e.g.
     *            seqan3::qualified, seqan3::masked, seqan3::structure_rna and
     *            seqan3::structure_aa).
     */
    //!\brief Construction from initialiser-list.
    constexpr alphabet_tuple_base(component_types... components) noexcept
    {
        assign_rank(rank_sum_helper(components..., std::make_index_sequence<sizeof...(component_types)>{}));
    }

    /*!\brief Construction via a value of one of the components.
     * \tparam component_type Must be one uniquely contained in the type list of the composite.
     * \param  alph           The value of a component that should be assigned.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_value_construction.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    template <typename component_type>
        requires (!std::is_base_of_v<alphabet_tuple_base, component_type>) && is_unique_component<component_type>
    constexpr explicit alphabet_tuple_base(component_type const alph) noexcept : alphabet_tuple_base{}
    {
        get<component_type>(*this) = alph;
    }

    /*!\brief Construction via a value of a subtype that is assignable to one of the components.
     * \tparam indirect_component_type Type that models seqan3::weakly_assignable_from for
     *                                 one of the component types.
     * \param  alph                    The value that should be assigned.
     *
     * Note that the value will be assigned to the **FIRST** type T that fulfils
     * `assignable_from<T, indirect_component_type>`, regardless if other types are also
     * fit for assignment.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_subtype_construction.cpp
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename indirect_component_type>
        requires ((detail::instantiate_if_v<
                       detail::lazy<std::is_convertible, indirect_component_type, component_types>,
                       detail::tuple_general_guard<derived_type, indirect_component_type, component_types...>>
                   || ...))
    constexpr explicit alphabet_tuple_base(indirect_component_type const alph) noexcept : alphabet_tuple_base{}
    {
        using component_predicate = detail::implicitly_convertible_from<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        component_type tmp(alph); // delegate construction
        get<component_type>(*this) = tmp;
    }

    //!\cond
    template <typename indirect_component_type>
        requires ((!(detail::instantiate_if_v<
                         detail::lazy<std::is_convertible, indirect_component_type, component_types>,
                         detail::tuple_general_guard<derived_type, indirect_component_type, component_types...>>
                     || ...))
                  && (detail::instantiate_if_v<
                          detail::lazy<std::is_constructible, component_types, indirect_component_type>,
                          detail::tuple_general_guard<derived_type, indirect_component_type, component_types...>>
                      || ...))
    constexpr explicit alphabet_tuple_base(indirect_component_type const alph) noexcept : alphabet_tuple_base{}
    {
        using component_predicate = detail::constructible_from<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        component_type tmp(alph); // delegate construction
        get<component_type>(*this) = tmp;
    }
    //!\endcond

    /*!\brief Assignment via a value of one of the components.
     * \tparam component_type One of the component types. Must be uniquely contained in the type list of the composite.
     * \param  alph           The value of a component that should be assigned.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_value_assignment.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    template <typename component_type>
        requires (!std::derived_from<component_type, alphabet_tuple_base>) && is_unique_component<component_type>
    constexpr derived_type & operator=(component_type const alph) noexcept
    {
        get<component_type>(*this) = alph;
        return static_cast<derived_type &>(*this);
    }

    /*!\brief Assignment via a value of a subtype that is assignable to one of the components.
     * \tparam indirect_component_type Type that models seqan3::weakly_assignable_from for
     *                                 one of the component types.
     * \param  alph                    The value of a component that should be assigned.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_subtype_assignment.cpp
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename indirect_component_type>
        requires ((!std::derived_from<indirect_component_type, alphabet_tuple_base>)
                  && (!is_unique_component<indirect_component_type>)
                  && (std::assignable_from<component_types, indirect_component_type> || ...))
    constexpr derived_type & operator=(indirect_component_type const alph) noexcept
    {
        using component_predicate = detail::assignable_from<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        get<component_type>(*this) = alph; // delegate assignment
        return static_cast<derived_type &>(*this);
    }
    //!\cond
    // If not assignable but implicit convertible, convert first and assign afterwards
    template <typename indirect_component_type>
        requires ((!std::derived_from<indirect_component_type, alphabet_tuple_base>)
                  && (!is_unique_component<indirect_component_type>)
                  && (!(std::assignable_from<component_types, indirect_component_type> || ...))
                  && (std::convertible_to<indirect_component_type, component_types> || ...))
    constexpr derived_type & operator=(indirect_component_type const alph) noexcept
    {
        using component_predicate = detail::implicitly_convertible_from<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        component_type tmp(alph);
        get<component_type>(*this) = tmp;
        return static_cast<derived_type &>(*this);
    }

    template <typename indirect_component_type>
        requires ((!std::derived_from<indirect_component_type, alphabet_tuple_base>)
                  && (!is_unique_component<indirect_component_type>)
                  && (!(std::assignable_from<component_types, indirect_component_type> || ...))
                  && (!(std::convertible_to<indirect_component_type, component_types> || ...))
                  && (std::constructible_from<component_types, indirect_component_type> || ...))
    constexpr derived_type & operator=(indirect_component_type const alph) noexcept
    {
        using component_predicate = detail::constructible_from<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        component_type tmp(alph); // delegate construction
        get<component_type>(*this) = tmp;
        return static_cast<derived_type &>(*this);
    }
    //!\endcond
    //!\}

    /*!\name Read functions
     * \brief All read operations are constant complexity.
     * \{
     */
    /*!\brief Tuple-like access to the contained components.
     * \tparam index Return the i-th element.
     * \returns A proxy to the contained element that models the same alphabets concepts and supports assignment.
     *
     * \stableapi{Since version 3.1.}
     */
    template <size_t index>
    friend constexpr auto get(alphabet_tuple_base & l) noexcept
    {
        static_assert(index < sizeof...(component_types), "Index out of range.");

        using t = seqan3::list_traits::at<index, component_list>;
        t val{};

        seqan3::assign_rank_to(l.to_component_rank<index>(), val);

        return component_proxy<t, index>{val, l};
    }

    /*!\copybrief seqan3::alphabet_tuple_base::get
     * \tparam type Return the element of specified type; only available if the type is unique in the set of components.
     * \returns A proxy to the contained element that models the same alphabets concepts and supports assignment.
     *
     * \stableapi{Since version 3.1.}
     */
    template <typename type>
    friend constexpr auto get(alphabet_tuple_base & l) noexcept
        requires is_unique_component<type>
    {
        return get<seqan3::list_traits::find<type, component_list>>(l);
    }

    /*!\copybrief seqan3::alphabet_tuple_base::get
     * \tparam index Return the i-th element.
     * \returns A copy of the contained element.
     *
     * \stableapi{Since version 3.1.}
     */
    template <size_t index>
    friend constexpr auto get(alphabet_tuple_base const & l) noexcept
    {
        static_assert(index < sizeof...(component_types), "Index out of range.");

        using t = seqan3::list_traits::at<index, component_list>;

        return seqan3::assign_rank_to(l.to_component_rank<index>(), t{});
    }

    /*!\copybrief seqan3::alphabet_tuple_base::get
     * \tparam type Return the element of specified type; only available if the type is unique in the set of components.
     * \returns A copy of the contained element.
     *
     * \stableapi{Since version 3.1.}
     */
    template <typename type>
    friend constexpr type get(alphabet_tuple_base const & l) noexcept
        requires is_unique_component<type>
    {
        return get<seqan3::list_traits::find<type, component_list>>(l);
    }

    /*!\brief Implicit cast to a single letter. Works only if the type is unique in the type list.
     *
     * \stableapi{Since version 3.1.}
     */
    template <typename type>
    constexpr operator type() const noexcept
        requires is_unique_component<type>
    {
        return get<type>(*this);
    }
    //!\}

    /*!\name Comparison operators (against indirect component_list)
     * \brief These overloads enable comparison for all types that a component type is comparable with.
     * \{
     */
    /*!\brief Comparison against types comparable with components.
     * \tparam indirect_component_type Must be comparable with a component's type.
     * \param lhs Left-hand-side of comparison.
     * \param rhs Right-hand-side of comparison.
     * \returns `true` or `false`.
     *
     * \details
     *
     * To determine (in-)equality/order, it is first deduced which component the argument is comparable with.
     * The tuple is then cast to that type and the resulting value compared with the argument.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator==(derived_type_t const lhs, indirect_component_type const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_eq_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        using component_predicate = detail::weakly_equality_comparable_with_<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        return get<component_type>(lhs) == rhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator==(indirect_component_type const lhs, derived_type_t const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_eq_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        return rhs == lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator!=(derived_type_t const lhs, indirect_component_type const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_eq_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        using component_predicate = detail::weakly_equality_comparable_with_<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        return get<component_type>(lhs) != rhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator!=(indirect_component_type const lhs, derived_type_t const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_eq_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        return rhs != lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator<(derived_type_t const lhs, indirect_component_type const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        using component_predicate = detail::weakly_ordered_with_<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        return get<component_type>(lhs) < rhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator<(indirect_component_type const lhs, derived_type_t const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        return rhs > lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator<=(derived_type_t const lhs, indirect_component_type const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        using component_predicate = detail::weakly_ordered_with_<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        return get<component_type>(lhs) <= rhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator<=(indirect_component_type const lhs, derived_type_t const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        return rhs >= lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator>(derived_type_t const lhs, indirect_component_type const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        using component_predicate = detail::weakly_ordered_with_<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        return get<component_type>(lhs) > rhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator>(indirect_component_type const lhs, derived_type_t const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        return rhs < lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator>=(derived_type_t const lhs, indirect_component_type const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        using component_predicate = detail::weakly_ordered_with_<indirect_component_type>;
        constexpr auto component_position =
            seqan3::list_traits::find_if<component_predicate::template invoke, component_list>;
        using component_type = seqan3::list_traits::at<component_position, component_list>;
        return get<component_type>(lhs) >= rhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::operator==(derived_type_t const lhs, indirect_component_type const rhs)
    template <typename derived_type_t, typename indirect_component_type>
    friend constexpr auto operator>=(indirect_component_type const lhs, derived_type_t const rhs) noexcept
        -> std::enable_if_t<
            detail::tuple_order_guard<derived_type_t, derived_type, indirect_component_type, component_types...>,
            bool>
    {
        return rhs <= lhs;
    }
    //!\}

private:
    //!\brief Return the rank of the i-th component.
    template <size_t index>
    constexpr rank_type to_component_rank() const noexcept
    {
        if constexpr (alphabet_size < 1024) // computation is cached for small alphabets
        {
            return rank_to_component_rank[index][to_rank()];
        }
        else
        {
            return (to_rank() / cummulative_alph_sizes[index])
                 % seqan3::alphabet_size<pack_traits::at<index, component_types...>>;
        }
    }

    //!\brief Assign via the rank of i-th component (does not update other components' state).
    template <size_t index>
    constexpr void assign_component_rank(ptrdiff_t const r) noexcept
    {
        assign_rank(static_cast<ptrdiff_t>(to_rank())
                    + ((r - static_cast<ptrdiff_t>(to_component_rank<index>()))
                       * static_cast<ptrdiff_t>(cummulative_alph_sizes[index])));
    }

    //!\brief For the given components, compute the combined rank.
    template <std::size_t... idx>
    static constexpr rank_type rank_sum_helper(component_types... components,
                                               std::index_sequence<idx...> const &) noexcept
    {
        return ((seqan3::to_rank(components) * cummulative_alph_sizes[idx]) + ...);
    }

    //!\brief The cumulative alphabet size products are cached.
    static constexpr std::array<rank_type, component_list::size()> cummulative_alph_sizes{
        []() constexpr
        {
            // create array (1, |sigma1|, |sigma1|*|sigma2|,  ... ,  |sigma1|*...*|sigmaN|)
            std::array<rank_type, component_list::size() + 1> ret{};
            ret[0] = 1;

            size_t count = 1;
            using reverse_list_t = decltype(seqan3::list_traits::detail::reverse(component_list{}));

            seqan3::detail::for_each<reverse_list_t>(
                [&](auto alphabet_type_identity) constexpr
                {
                    using alphabet_t = typename decltype(alphabet_type_identity)::type;
                    ret[count] = static_cast<rank_type>(seqan3::alphabet_size<alphabet_t> * ret[count - 1]);
                    ++count;
                });

            // reverse and strip one: (|sigma1|*...*|sigmaN-1|, ... |sigma1|*|sigma2|, |sigma1|, 1)
            // reverse order guarantees that the first alphabet is the most significant rank contributer
            // resulting in element-wise alphabetical ordering on comparison
            std::array<rank_type, component_list::size()> ret2{};

            for (size_t i = 0; i < component_list::size(); ++i)
                ret2[i] = ret[component_list::size() - i - 1];

            return ret2;
        }()};

    //!\brief Conversion table from rank to the i-th component's rank.
    static constexpr std::array < std::array<rank_type,
                                             alphabet_size<1024u ? alphabet_size : 0u>, // not for big alphs
                                             list_traits::size<component_list>>
            rank_to_component_rank{
                []() constexpr
                {
                    std::array < std::array<rank_type,
                                            alphabet_size<1024u ? alphabet_size : 0u>, // not for big alphs
                                            list_traits::size<component_list>>
                            ret{};

                    if constexpr (alphabet_size < 1024u)
                    {
                        std::array<size_t, alphabet_size> alph_sizes{seqan3::alphabet_size<component_types>...};

                        for (size_t i = 0; i < list_traits::size<component_list>; ++i)
                            for (size_t j = 0; j < static_cast<size_t>(alphabet_size); ++j)
                                ret[i][j] = (j / cummulative_alph_sizes[i]) % alph_sizes[i];
                    }

                    return ret;
                }()};
};

/*!\brief Specialisation of seqan3::alphabet_proxy that updates the rank of the alphabet_tuple_base.
 * \tparam alphabet_type The type of the emulated component.
 * \tparam index         The index of the emulated component.
 *
 * \noapi
 */
template <typename derived_type, typename... component_types>
    requires (detail::writable_constexpr_semialphabet<component_types> && ...) && (std::regular<component_types> && ...)
template <typename alphabet_type, size_t index>
class alphabet_tuple_base<derived_type, component_types...>::component_proxy :
    public alphabet_proxy<component_proxy<alphabet_type, index>, alphabet_type>
{
private:
    //!\brief The base type.
    using base_t = alphabet_proxy<component_proxy<alphabet_type, index>, alphabet_type>;
    //!\brief Befriend the base type so it can call our seqan3::alphabet_tuple_base::component_proxy::on_update().
    friend base_t;

    //!\brief Store a pointer to the parent object so we can update it.
    alphabet_tuple_base * parent;

    //!\brief The implementation updates the rank in the parent object.
    constexpr void on_update() noexcept
    {
        parent->assign_component_rank<index>(this->to_rank());
    }

public:
    //Import from base type:
    using base_t::operator=;

    /*!\name Constructors, destructor and assignment
        * \{
        */
    //!\brief Deleted, because using this proxy without parent would be undefined behaviour.
    component_proxy() = delete;
    constexpr component_proxy(component_proxy const &) = default;             //!< Defaulted.
    constexpr component_proxy(component_proxy &&) = default;                  //!< Defaulted.
    constexpr component_proxy & operator=(component_proxy const &) = default; //!< Defaulted.
    constexpr component_proxy & operator=(component_proxy &&) = default;      //!< Defaulted.
    ~component_proxy() = default;                                             //!< Defaulted.

    //!\brief Construct from an alphabet letter and reference to the parent object.
    constexpr component_proxy(alphabet_type const l, alphabet_tuple_base & p) : base_t{l}, parent{&p}
    {}

    // Does not inherit the base's constructor for alphabet_type so as not to cause ambiguity
    //!\}

    /*!\name Comparison operators (proxy type against parent)
        * \brief Comparison against the seqan3::alphabet_tuple_base that this proxy originates from (necessary
        *        to prevent recursive template instantiation in the tuple).
        * \{
        */
    /*!\brief Comparison against the alphabet tuple by casting self and tuple to the emulated type.
        * \param lhs Left-hand-side of comparison.
        * \param rhs Right-hand-side of comparison.
        */
    friend constexpr bool operator==(derived_type const lhs, component_proxy const rhs) noexcept
    {
        return get<index>(lhs) == static_cast<alphabet_type>(rhs);
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator==(component_proxy<alphabet_type, index> const lhs, derived_type const rhs) noexcept
    {
        return rhs == lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator!=(derived_type const lhs, component_proxy const rhs) noexcept
    {
        return get<index>(lhs) != static_cast<alphabet_type>(rhs);
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator!=(component_proxy<alphabet_type, index> const lhs, derived_type const rhs) noexcept
    {
        return rhs != lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator<(derived_type const lhs, component_proxy const rhs) noexcept
    {
        return get<index>(lhs) < static_cast<alphabet_type>(rhs);
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator<(component_proxy<alphabet_type, index> const lhs, derived_type const rhs) noexcept
    {
        return rhs > lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator<=(derived_type const lhs, component_proxy const rhs) noexcept
    {
        return get<index>(lhs) <= static_cast<alphabet_type>(rhs);
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator<=(component_proxy<alphabet_type, index> const lhs, derived_type const rhs) noexcept
    {
        return rhs >= lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator>(derived_type const lhs, component_proxy const rhs) noexcept
    {
        return get<index>(lhs) > static_cast<alphabet_type>(rhs);
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator>(component_proxy<alphabet_type, index> const lhs, derived_type const rhs) noexcept
    {
        return rhs < lhs;
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator>=(derived_type const lhs, component_proxy const rhs) noexcept
    {
        return get<index>(lhs) >= static_cast<alphabet_type>(rhs);
    }

    //!\copydoc seqan3::alphabet_tuple_base::component_proxy::operator==(derived_type const, component_proxy const)
    friend constexpr bool operator>=(component_proxy<alphabet_type, index> const lhs, derived_type const rhs) noexcept
    {
        return rhs <= lhs;
    }
    //!\}
};

} // namespace seqan3

namespace std
{

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::transformation_trait
 * \ingroup alphabet_composite
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 *
 * \stableapi{Since version 3.1.}
 */
template <std::size_t i, seqan3::detail::alphabet_tuple_like tuple_t>
struct tuple_element<i, tuple_t>
{
    //!\brief Element type.
    using type = seqan3::list_traits::at<i, typename tuple_t::seqan3_required_types>;
};

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::unary_type_trait
 * \ingroup alphabet_composite
 * \see std::tuple_size_v
 *
 * \stableapi{Since version 3.1.}
 */
template <seqan3::detail::alphabet_tuple_like tuple_t>
struct tuple_size<tuple_t> :
    public std::integral_constant<size_t, seqan3::list_traits::size<typename tuple_t::seqan3_required_types>>
{};

} // namespace std
