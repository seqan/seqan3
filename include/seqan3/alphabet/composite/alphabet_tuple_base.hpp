// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_tuple_base.
 */

#pragma once

#include <cassert>
#include <utility>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/composite/detail.hpp>
#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/detail/alphabet_proxy.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/core/type_traits/pack.hpp>
#include <seqan3/core/type_traits/transformation_trait_or.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief Evaluates to true if one of the components of seqan3::alphabet_tuple_base satisifes a compile-time
 *        predicate.
 * \tparam tuple_t         A specialisation of seqan3::alphabet_tuple_base.
 * \tparam tuple_derived_t The CRTP derived type of `tuple_t`.
 * \tparam fun_t               A template template that takes target_t as argument and exposes an `invoke` member type
 *                             that evaluates some predicate and returns `std::true_type` or `std::false_type`.
 * \tparam target_t            The type you wish to query.
 * \ingroup composite
 *
 * \details
 *
 * To prevent recursive template and/or concept instantiation this call needs to be guarded against many exceptions.
 * See the source file for more details.
 */

// anchor is false
template <typename tuple_t, typename tuple_derived_t, template <typename> typename fun_t, typename other_t>
inline bool constexpr one_component_is = false;

//!\cond

// default
template <typename ... tuple_comps,
          typename tuple_derived_t,
          template <typename> typename fun_t,
          typename other_t>
inline bool constexpr one_component_is<alphabet_tuple_base<tuple_derived_t, tuple_comps...>,
                                       tuple_derived_t,
                                       fun_t,
                                       other_t> =
    !meta::empty<meta::find_if<meta::list<tuple_comps...>, fun_t<other_t>>>::value;
    //TODO do without meta

// guard against self
template <typename ... tuple_comps,
          typename tuple_derived_t,
          template <typename> typename fun_t>
inline bool constexpr one_component_is<alphabet_tuple_base<tuple_derived_t, tuple_comps...>,
                                       tuple_derived_t,
                                       fun_t,
                                       alphabet_tuple_base<tuple_derived_t, tuple_comps...>> = false;

// guard against self (derived)
template <typename ... tuple_comps,
          typename tuple_derived_t,
          template <typename> typename fun_t>
inline bool constexpr one_component_is<alphabet_tuple_base<tuple_derived_t, tuple_comps...>,
                                       tuple_derived_t,
                                       fun_t,
                                       tuple_derived_t> = false;

// guard against types that have conversion operators to derived
template <typename ... tuple_comps,
          typename tuple_derived_t,
          template <typename> typename fun_t,
          typename other_t>
    requires ConvertibleToByMember<other_t, tuple_derived_t>
inline bool constexpr one_component_is<alphabet_tuple_base<tuple_derived_t, tuple_comps...>,
                                       tuple_derived_t,
                                       fun_t,
                                       other_t> = false;

// guard against components
template <typename ... tuple_comps,
          typename tuple_derived_t,
          template <typename> typename fun_t,
          typename other_t>
    requires type_in_pack_v<other_t, tuple_comps...>
//     requires meta::in<recursive_tuple_components<alphabet_tuple_base<tuple_derived_t, tuple_comps...>>::type,
//                       other_t>::value
inline bool constexpr one_component_is<alphabet_tuple_base<tuple_derived_t, tuple_comps...>,
                                       tuple_derived_t,
                                       fun_t,
                                       other_t> = false;

// during comparisons, guard against types that could be converted to self (because that is preferred)
// (may not be done during assignment or construction because of recursiveness)
template <typename ... tuple_comps,
          typename tuple_derived_t,
          typename other_t>
    requires ImplicitlyConvertibleTo<other_t, tuple_derived_t>
inline bool constexpr one_component_is<alphabet_tuple_base<tuple_derived_t, tuple_comps...>,
                                       tuple_derived_t,
                                       weakly_equality_comparable_with,
                                       other_t> = false;
template <typename ... tuple_comps,
          typename tuple_derived_t,
          typename other_t>
    requires ImplicitlyConvertibleTo<other_t, tuple_derived_t>
inline bool constexpr one_component_is<alphabet_tuple_base<tuple_derived_t, tuple_comps...>,
                                       tuple_derived_t,
                                       weakly_ordered_with,
                                       other_t> = false;
//!\endcond

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
 * \ingroup composite
 * \implements seqan3::WritableSemialphabet
 * \if DEV
 * \implements seqan3::detail::ConstexprWritableSemialphabet
 * \tparam component_types Types of letters; must model seqan3::detail::WritableConstexprSemialphabet.
 * \else
 * \tparam component_types Types of letters; must model seqan3::WritableSemialphabet and all required function calls
 * need to be callable in `constexpr`-context.
 * \endif
 *
 *
 * This data structure is a CRTP base class for combined alphabets, where the different
 * alphabet letters exist independently as a components, similar to a tuple.
 *
 * Short description:
 *   * combines multiple alphabets as independent components, similar to a tuple;
 *   * models seqan3::TupleLike, i.e. provides a get interface to its component_list;
 *   * is itself a seqan3::WritableSemialphabet, but most derived types implement the full seqan3::WritableAlphabet;
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
 */
template <typename derived_type,
          typename ...component_types>
//!\cond
    requires (detail::WritableConstexprSemialphabet<component_types> && ...) &&
             (!std::is_reference_v<component_types> && ...)
//!\endcond
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

    //!\brief A meta::list The types of each component in the composite
    using component_list = meta::list<component_types...>;

    //!\brief Is set to `true` if the type is contained in the type list.
    template <typename type>
    static constexpr bool is_component =
        meta::in<component_list, type>::value;

    //!\brief Is set to `true` if the type is uniquely contained in the type list.
    template <typename type>
    static constexpr bool is_unique_component =
        is_component<type> &&
        (meta::find_index<component_list, type>::value == meta::reverse_find_index<component_list, type>::value);

    /*!\brief Specialisation of seqan3::alphabet_proxy that updates the rank of the alphabet_tuple_base.
     * \tparam alphabet_type The type of the emulated component.
     * \tparam index         The index of the emulated component.
     */
    template <typename alphabet_type, size_t index>
    class component_proxy : public alphabet_proxy<component_proxy<alphabet_type, index>, alphabet_type>
    {
    private:
        //!\brief The base type.
        using base_t = alphabet_proxy<component_proxy<alphabet_type, index>, alphabet_type>;
        //!\brief Befriend the base type so it can call our #on_update().
        friend base_t;

        //!\brief Store a pointer to the parent object so we can update it.
        alphabet_tuple_base *parent;

        //!\brief The implementation updates the rank in the parent object.
        constexpr void on_update() noexcept
        {
            parent->assign_rank(
                parent->to_rank()
                - parent->template to_component_rank<index>() * alphabet_tuple_base::cummulative_alph_sizes[index]
                + to_rank() * alphabet_tuple_base::cummulative_alph_sizes[index]);
        }

        /*!\name Associated types
         * \{
         */
        using typename base_t::rank_type;
        using typename base_t::char_type;
        using typename base_t::phred_type;
        //!\}

    public:
        //Import from base type:
        using base_t::to_rank;
        using base_t::alphabet_size;
        using base_t::operator=;

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr component_proxy() : base_t{}, parent{} {}                        //!< Defaulted.
        constexpr component_proxy(component_proxy const &) = default;              //!< Defaulted.
        constexpr component_proxy(component_proxy &&) = default;                   //!< Defaulted.
        constexpr component_proxy & operator=(component_proxy const &) = default;  //!< Defaulted.
        constexpr component_proxy & operator=(component_proxy &&) = default;       //!< Defaulted.
        ~component_proxy() = default;                                              //!< Defaulted.

        constexpr component_proxy(alphabet_type const l, alphabet_tuple_base & p) :
            base_t{l}, parent{&p}
        {}

        // Does not inherit the base's constructor for alphabet_type so as not to cause ambiguity
        //!\}
    };

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_tuple_base() noexcept : base_t{} {}
    constexpr alphabet_tuple_base(alphabet_tuple_base const &) = default;
    constexpr alphabet_tuple_base(alphabet_tuple_base &&) = default;
    constexpr alphabet_tuple_base & operator=(alphabet_tuple_base const &) = default;
    constexpr alphabet_tuple_base & operator=(alphabet_tuple_base &&) = default;
    ~alphabet_tuple_base() = default;

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
    using base_t::to_rank;
    using base_t::assign_rank;

    //!\brief Export this type's components in a visible manner.
    //!\private
    using seqan3_tuple_components = component_list;
    //!\brief Export this type's components and possibly the components' components in a visible manner.
    //!\private
    using seqan3_recursive_tuple_components =
        meta::concat<component_list,
                     detail::transformation_trait_or_t<detail::recursive_tuple_components<component_types>,
                                                       meta::list<>>...>;

    /*!\name Constructors, destructor and assignment
     * \{
     * \attention Please do not directly use the CRTP base class. The functions
     *            are only public for the usage in their derived classes (e.g.
     *            seqan3::qualified, seqan3::masked, seqan3::structure_rna and
     *            seqan3::structure_aa).
     */
    //!\brief Construction from initialiser-list.
    constexpr alphabet_tuple_base(component_types ... components) noexcept
    {
        assign_rank(rank_sum_helper(components..., std::make_index_sequence<sizeof...(component_types)>{}));
    }

    /*!\brief Construction via a value of one of the components.
     * \tparam component_type Must be one uniquely contained in the type
                              list of the composite.
     * \param  alph           The value of a component that should be assigned.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_value_construction.cpp
     */
    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr explicit alphabet_tuple_base(component_type const alph) noexcept : alphabet_tuple_base{}
    {
        get<component_type>(*this) = alph;
    }

    /*!\brief Construction via a value of a subtype that is assignable to one of the components.
     * \tparam indirect_component_type Type that models seqan3::WeaklyAssignable for
     *                                 one of the component types.
     * \param  alph                    The value that should be assigned.
     *
     * Note that the value will be assigned to the **FIRST** type T that fulfils
     * `Assignable<T, indirect_component_type>`, regardless if other types are also
     * fit for assignment.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_subtype_construction.cpp
     */
    template <typename indirect_component_type>
    //!\cond
       requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::implicitly_convertible_from, indirect_component_type>
    //!\endcond
    constexpr explicit alphabet_tuple_base(indirect_component_type const alph) noexcept : alphabet_tuple_base{}
    {
       using component_type = meta::front<meta::find_if<component_list, detail::implicitly_convertible_from<indirect_component_type>>>;
       component_type tmp(alph); // delegate construction
       get<component_type>(*this) = tmp;
    }

    //!\cond
    template <typename indirect_component_type>
       requires !detail::one_component_is<alphabet_tuple_base, derived_type, detail::implicitly_convertible_from, indirect_component_type> &&
                 detail::one_component_is<alphabet_tuple_base, derived_type, detail::constructible_from, indirect_component_type>
    constexpr explicit alphabet_tuple_base(indirect_component_type const alph) noexcept : alphabet_tuple_base{}
    {
       using component_type = meta::front<meta::find_if<component_list, detail::constructible_from<indirect_component_type>>>;
       component_type tmp(alph); // delegate construction
       get<component_type>(*this) = tmp;
    }
    //!\endcond

    /*!\brief Assignment via a value of one of the components.
     * \tparam component_type One of the component types. Must be uniquely
     *                        contained in the type list of the composite.
     * \param  alph           The value of a component that should be assigned.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_value_assignment.cpp
     */
    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr derived_type & operator=(component_type const alph) noexcept
    {
        get<component_type>(*this) = alph;
        return static_cast<derived_type &>(*this);
    }

    /*!\brief Assignment via a value of a subtype that is assignable to one of the components.
     * \tparam indirect_component_type Type that models seqan3::WeaklyAssignable for
     *                                 one of the component types.
     * \param  alph                    The value of a component that should be assigned.
     *
     * Note: Since the alphabet_tuple_base is a CRTP base class, we show the working examples
     * with one of its derived classes (seqan3::qualified).
     * \include test/snippet/alphabet/composite/alphabet_tuple_base_subtype_assignment.cpp
     */
    template <typename indirect_component_type>
    //!\cond
        requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::assignable_from, indirect_component_type>
    //!\endcond
    constexpr derived_type & operator=(indirect_component_type const alph) noexcept
    {
        using component_type = meta::front<meta::find_if<component_list, detail::assignable_from<indirect_component_type>>>;
        get<component_type>(*this) = alph; // delegate assignment
        return static_cast<derived_type &>(*this);
    }
    //!\cond
    // If not assignable but implicit convertible, convert first and assign afterwards
    template <typename indirect_component_type>
    //!\cond
        requires !detail::one_component_is<alphabet_tuple_base, derived_type, detail::assignable_from, indirect_component_type> &&
                 detail::one_component_is<alphabet_tuple_base, derived_type, detail::implicitly_convertible_from, indirect_component_type>
    //!\endcond
    constexpr derived_type & operator=(indirect_component_type const alph) noexcept
    {
        using component_type = meta::front<meta::find_if<component_list, detail::implicitly_convertible_from<indirect_component_type>>>;
        component_type tmp(alph);
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
     */
    template <size_t index>
    friend constexpr auto get(alphabet_tuple_base & l) noexcept
    {
        static_assert(index < sizeof...(component_types), "Index out of range.");

        using t = meta::at_c<component_list, index>;
        t val{};

        using seqan3::assign_rank_to;
        assign_rank_to(l.to_component_rank<index>(), val);

        return component_proxy<t, index>{val, l};
    }

    /*!\copybrief get
     * \tparam type Return the element of specified type; only available if the type is unique in the set of components.
     * \returns A proxy to the contained element that models the same alphabets concepts and supports assignment.
     */
    template <typename type>
    friend constexpr auto get(alphabet_tuple_base & l) noexcept
    //!\cond
        requires is_unique_component<type>
    //!\endcond
    {
        return get<meta::find_index<component_list, type>::value>(l);
    }

    /*!\copybrief get
     * \tparam index Return the i-th element.
     * \returns A copy of the contained element.
     */
    template <size_t index>
    friend constexpr auto get(alphabet_tuple_base const & l) noexcept
    {
        static_assert(index < sizeof...(component_types), "Index out of range.");

        using t = meta::at_c<component_list, index>;
        t val{};

        using seqan3::assign_rank_to;
        assign_rank_to(l.to_component_rank<index>(), val);

        return val;
    }

    /*!\copybrief get
     * \tparam type Return the element of specified type; only available if the type is unique in the set of components.
     * \returns A copy of the contained element.
     */
    template <typename type>
    friend constexpr type get(alphabet_tuple_base const & l) noexcept
    //!\cond
        requires is_unique_component<type>
    //!\endcond
    {
        return get<meta::find_index<component_list, type>::value>(l);
    }

    /*!\brief Implicit cast to a single letter. Works only if the type is unique in the type list.
     */
    template <typename type>
    constexpr operator type() const noexcept
    //!\cond
        requires is_unique_component<type>
    //!\endcond
    {
        return get<type>(*this);
    }
    //!\}

    /*!\name Comparison operators (against indirect component_list)
     * \brief `*this` is cast to the component type before comparison. (These overloads enable comparison for all types
     *        that a component type is comparable with).
     * \{
     */
    template <typename indirect_component_type>
    constexpr bool operator==(indirect_component_type const rhs) const noexcept
    //!\cond
        requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::weakly_equality_comparable_with, indirect_component_type>
    //!\endcond
    {
        using component_type = meta::front<meta::find_if<component_list, detail::weakly_equality_comparable_with<indirect_component_type>>>;
        return get<component_type>(*this) == rhs;
    }

    template <typename indirect_component_type>
    constexpr bool operator!=(indirect_component_type const rhs) const noexcept
    //!\cond
        requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::weakly_equality_comparable_with, indirect_component_type>
    //!\endcond
    {
        using component_type = meta::front<meta::find_if<component_list, detail::weakly_equality_comparable_with<indirect_component_type>>>;
        return get<component_type>(*this) != rhs;
    }

    template <typename indirect_component_type>
    constexpr bool operator<(indirect_component_type const rhs) const noexcept
    //!\cond
        requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::weakly_ordered_with, indirect_component_type>
    //!\endcond
    {
        using component_type = meta::front<meta::find_if<component_list, detail::weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) < rhs;
    }

    template <typename indirect_component_type>
    constexpr bool operator>(indirect_component_type const rhs) const noexcept
    //!\cond
        requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::weakly_ordered_with, indirect_component_type>
    //!\endcond
    {
        using component_type = meta::front<meta::find_if<component_list, detail::weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) > rhs;
    }

    template <typename indirect_component_type>
    constexpr bool operator<=(indirect_component_type const rhs) const noexcept
    //!\cond
        requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::weakly_ordered_with, indirect_component_type>
    //!\endcond
    {
        using component_type = meta::front<meta::find_if<component_list, detail::weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) <= rhs;
    }

    template <typename indirect_component_type>
    constexpr bool operator>=(indirect_component_type const rhs) const noexcept
    //!\cond
        requires detail::one_component_is<alphabet_tuple_base, derived_type, detail::weakly_ordered_with, indirect_component_type>
    //!\endcond
    {
        using component_type = meta::front<meta::find_if<component_list, detail::weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) >= rhs;
    }
    //!\}

private:
    //!\brief Return the rank of the i-th component.
    template <size_t index>
    constexpr rank_type to_component_rank() const noexcept
    {
        return (to_rank() / cummulative_alph_sizes[index]) % seqan3::alphabet_size<meta::at_c<component_list, index>>;
    }

    //!\brief The cumulative alphabet size products are cached.
    static constexpr std::array<rank_type, component_list::size()> cummulative_alph_sizes
    {
        [] () constexpr
        {
            // create array (1, |sigma1|, |sigma1|*|sigma2|,  ... ,  |sigma1|*...*|sigmaN|)
            std::array<rank_type, component_list::size() + 1> ret{};
            ret[0] = 1;
            size_t count = 1;
            meta::for_each(meta::reverse<component_list>{}, [&] (auto alph) constexpr
            {
                ret[count] = static_cast<rank_type>(seqan3::alphabet_size<decltype(alph)> * ret[count - 1]);
                ++count;
            });

            // reverse and strip one: (|sigma1|*...*|sigmaN-1|, ... |sigma1|*|sigma2|, |sigma1|, 1)
            // reverse order guarantees that the first alphabet is the most significant rank contributer
            // resulting in element-wise alphabetical ordering on comparison
            std::array<rank_type, component_list::size()> ret2{};
            for (size_t i = 0; i < component_list::size(); ++i)
                ret2[i] = ret[component_list::size() - i - 1];

            return ret2;
        }()
    };

    //!\brief For the given components, compute the combined rank.
    template <std::size_t ...idx>
    static constexpr rank_type rank_sum_helper(component_types ... components, std::index_sequence<idx...> const &) noexcept
    {
        using seqan3::to_rank;
        return ((to_rank(components) * cummulative_alph_sizes[idx]) + ...);
    }
};

/*!\name Comparison operators
 * \relates seqan3::alphabet_tuple_base
 * \{
 * \brief Free function comparison operators that forward to member operators (for types != self).
 */
template <typename indirect_component_type, typename derived_type, typename ...component_types>
//!\cond
    requires detail::WeaklyEqualityComparableByMembersWith<derived_type, indirect_component_type> &&
             !detail::WeaklyEqualityComparableByMembersWith<indirect_component_type, derived_type>
//!\endcond
constexpr bool operator==(indirect_component_type const lhs,
                          alphabet_tuple_base<derived_type, component_types...> const rhs) noexcept
{
    return rhs == lhs;
}

template <typename indirect_component_type, typename derived_type, typename ...indirect_component_types>
//!\cond
    requires detail::WeaklyEqualityComparableByMembersWith<derived_type, indirect_component_type> &&
             !detail::WeaklyEqualityComparableByMembersWith<indirect_component_type, derived_type>
//!\endcond
constexpr bool operator!=(indirect_component_type const lhs,
                          alphabet_tuple_base<derived_type, indirect_component_types...> const rhs) noexcept
{
    return rhs != lhs;
}

template <typename indirect_component_type, typename derived_type, typename ...indirect_component_types>
//!\cond
    requires detail::WeaklyOrderedByMembersWith<derived_type, indirect_component_type> &&
             !detail::WeaklyOrderedByMembersWith<indirect_component_type, derived_type>
//!\endcond
constexpr bool operator<(indirect_component_type const lhs,
                         alphabet_tuple_base<derived_type, indirect_component_types...> const rhs) noexcept
{
    return rhs > lhs;
}

template <typename indirect_component_type, typename derived_type, typename ...indirect_component_types>
//!\cond
    requires detail::WeaklyOrderedByMembersWith<derived_type, indirect_component_type> &&
             !detail::WeaklyOrderedByMembersWith<indirect_component_type, derived_type>
//!\endcond
constexpr bool operator>(indirect_component_type const lhs,
                         alphabet_tuple_base<derived_type, indirect_component_types...> const rhs) noexcept
{
    return rhs < lhs;
}

template <typename indirect_component_type, typename derived_type, typename ...indirect_component_types>
//!\cond
    requires detail::WeaklyOrderedByMembersWith<derived_type, indirect_component_type> &&
             !detail::WeaklyOrderedByMembersWith<indirect_component_type, derived_type>
//!\endcond
constexpr bool operator<=(indirect_component_type const lhs,
                          alphabet_tuple_base<derived_type, indirect_component_types...> const rhs) noexcept
{
    return rhs >= lhs;
}

template <typename indirect_component_type, typename derived_type, typename ...indirect_component_types>
//!\cond
    requires detail::WeaklyOrderedByMembersWith<derived_type, indirect_component_type> &&
             !detail::WeaklyOrderedByMembersWith<indirect_component_type, derived_type>
//!\endcond
constexpr bool operator>=(indirect_component_type const lhs,
                          alphabet_tuple_base<derived_type, indirect_component_types...> const rhs) noexcept
{
    return rhs <= lhs;
}
//!\}

} // namespace seqan3

namespace std
{

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::TransformationTrait
 * \ingroup composite
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <std::size_t i, seqan3::detail::AlphabetTupleBase tuple_t>
struct tuple_element<i, tuple_t>
{
    //!\brief Element type.
    using type = meta::at_c<typename tuple_t::seqan3_tuple_components, i>;
};

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::UnaryTypeTrait
 * \ingroup composite
 * \see std::tuple_size_v
 */
template <seqan3::detail::AlphabetTupleBase tuple_t>
struct tuple_size<tuple_t> :
    public std::integral_constant<size_t, tuple_t::seqan3_tuple_components::size()>
{};

} // namespace std
