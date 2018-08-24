// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains cartesian_composition.
 */

#pragma once

#include <cassert>
#include <utility>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief The CRTP base of alphabets that contain multiple (different) letters at one position.
 * \ingroup composition
 * \implements seqan3::semi_alphabet_concept
 * \tparam first_component_type Type of the first letter; must model seqan3::semi_alphabet_concept.
 * \tparam component_types      Types of further letters (up to 4); must model seqan3::semi_alphabet_concept.
 *
 * This data structure is CRTP base class for combined alphabets, where the different
 * alphabet letters exist independently, similar to a tuple. In fact this class
 * provides a tuple-like interface with `get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * \attention
 * This is a "pure base class", you cannot instantiate it, you can only inherit from it.
 * Most likely you are interested in using one of it's descendants like quality_composition.
 * \cond DEV
 * To make a derived class "complete", you should add at least the following:
 *   * .to_char() member
 *   * .assign_char() member
 *   * .operator=() members for all element types
 *   * a type deduction guide
 * \endcond
 *
 * \sa qualified
 * \sa masked
 */
template <typename derived_type,
          typename first_component_type,
          typename ...component_types>
//!\cond
    requires semi_alphabet_concept<first_component_type> && (semi_alphabet_concept<component_types> && ...)
//!\endcond
class cartesian_composition : public pod_tuple<first_component_type, component_types...>
{
private:
    //!\brief The base type of this class.
    using base_t = pod_tuple<first_component_type, component_types...>;

    //!\brief A meta::list The types of each component in the composition
    using components = meta::list<first_component_type, component_types...>;

    //!\brief Is set to `true` if the type is contained in the type list.
    template <typename type>
    static constexpr bool is_component =
        meta::in<components, type>::value;

    //!\brief Is set to `true` if the type is uniquely contained in the type list.
    template <typename type>
    static constexpr bool is_unique_component =
        is_component<type> &&
        (meta::find_index<components, type>::value == meta::reverse_find_index<components, type>::value);

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
      * Returns an std::true_type if the `type` is constructable from `T`.
      */
    template <typename T>
    struct constructible_from
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, std::is_constructible_v<type, T>>;
    };

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
      * Returns an std::true_type if the `T` is implicitly convertible to `type`.
      */
    template <typename T>
    struct implicitely_convertible_from
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, implicitly_convertible_to_concept<T, type>>;
    };

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
      * Returns an std::true_type if the `type` is assignable from `T`.
      */
    template <typename T>
    struct assignable_from
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, std::Assignable<type, T>>;
    };

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
      * Returns an std::true_type if the `type` is weakly equality comparable to `T`.
      */
    template <typename T>
    struct weakly_equality_comparable_with
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, std::detail::WeaklyEqualityComparableWith<type, T>>;
    };

    /*!\brief 'Callable' helper class that is invokable by meta::invoke.
      * Returns an std::true_type if the `type` is comparable via <,<=,>,>= to `T`.
      */
    template <typename T>
    struct weakly_ordered_with
    {
        //!\brief The returned type when invoked.
        template <typename type>
        using invoke = std::integral_constant<bool, weakly_ordered_with_concept<type, T>>;
    };

    //!\brief Is set to `true` if the one component type fulfils the FUN callable with `subtype`.
    template <template <typename> typename FUN, typename subtype>
    static constexpr bool one_component_is =
         !std::is_same_v<subtype, cartesian_composition> &&
         !std::is_same_v<subtype, derived_type> &&
         !is_component<subtype> &&
         !meta::empty<meta::find_if<components, FUN<subtype>>>::value;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    cartesian_composition() = default;
    constexpr cartesian_composition(cartesian_composition const &) = default;
    constexpr cartesian_composition(cartesian_composition &&) = default;
    constexpr cartesian_composition & operator=(cartesian_composition const &) = default;
    constexpr cartesian_composition & operator=(cartesian_composition &&) = default;
    ~cartesian_composition() = default;
    //!\}

public:
    //!\brief The type of value_size and `alphabet_size_v<cartesian_composition<...>>`
    using rank_type = detail::min_viable_uint_t<(alphabet_size_v<first_component_type> * ... *
                                                     alphabet_size_v<component_types>)>;

    //!\brief The product of the sizes of the individual alphabets.
    static constexpr rank_type value_size{(alphabet_size_v<first_component_type> * ... *
                                               alphabet_size_v<component_types>)};

    /*!\name Constructors, destructor and assignment
     * \{
     * \attention Please do not directly use the CRTP base class. The functions
     *            are only public for the usage in their derived classes (e.g.
     *            seqan3::qualified, seqan3::masked, seqan3::structure_rna and
     *            seqan3::structure_aa).
     */
    //!\brief Construction from initialiser-list.
    constexpr cartesian_composition(first_component_type first, component_types ... others) :
        base_t{first, others...}
    {}

    /*!\brief Construction via a value of one of the components.
     * \tparam component_type Must be one uniquely contained in the type
                              list of the composition.
     * \param  alph           The value of a component that should be assigned.
     *
     * ```cpp
     *     cartesian_composition<dna4, aa27> letter1{dna4::C}; // creates {dna4::C, aa27::A}
     *     cartesian_composition<dna4, aa27> letter2{aa27::K}; // creates {dna4::A, aa27::K}
     * ```
     */
    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr cartesian_composition(component_type const alph) : cartesian_composition{}
    {
        get<component_type>(*this) = alph;
    }

    /*!\brief Construction via a value of a subtype that is assignable to one of the components.
     * \tparam indirect_component_type Type that models the seqan3::is_assignable_concept for
     *                                 one of the component types.
     * \param  alph                    The value that should be assigned.
     *
     * Note that the value will be assigned to the **FIRST** type T that fulfils
     * the `assignable_concept<T, indirect_component_type>`, regardless if other types are also
     * fit for assignment.
     *
     * ```cpp
     *     // The following creates {dna4::C, aa27::A}
     *     cartesian_composition<gapped<dna4>, aa27> letter1{dna4::C};
     *     // The following creates {dna4::A, dna4::C}, since dna4 is also a type in the list
     *     cartesian_composition<gapped<dna4>, dna4> letter2{dna4::C};
     *     // The following creates {dna5::C, dna15::A}, since dna5 is the first type assignable to dna4
     *     cartesian_composition<dna5, dna15> letter2{dna4::C};
     * ```
     */
    template <typename indirect_component_type>
    //!\cond
       requires one_component_is<implicitely_convertible_from, indirect_component_type>
    //!\endcond
    constexpr cartesian_composition(indirect_component_type const alph) : cartesian_composition{}
    {
       using component_type = meta::front<meta::find_if<components, implicitely_convertible_from<indirect_component_type>>>;
       component_type tmp(alph); // delegate construction
       get<component_type>(*this) = tmp;
    }

    //!\cond
    template <typename indirect_component_type>
       requires !one_component_is<implicitely_convertible_from, indirect_component_type> &&
                 one_component_is<constructible_from, indirect_component_type>
    constexpr cartesian_composition(indirect_component_type const alph) : cartesian_composition{}
    {
       using component_type = meta::front<meta::find_if<components, constructible_from<indirect_component_type>>>;
       component_type tmp(alph); // delegate construction
       get<component_type>(*this) = tmp;
    }
    //!\endcond

    /*!\brief Assignment via a value of one of the components.
     * \tparam component_type One of the component types. Must be uniquely
     *                        contained in the type list of the composition.
     * \param  alph           The value of a component that should be assigned.
     *
     * ```cpp
     *     cartesian_composition<dna4, aa27> letter1{dna4::T, aa27::K};
     *
     *     letter1 = dna4::C // yields {dna4::C, aa27::K}
     *     letter1 = aa27::F // yields {dna4::C, aa27::F}
     * ```
     */
    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr derived_type & operator=(component_type const alph)
    {
        get<component_type>(*this) = alph;
        return static_cast<derived_type &>(*this);
    }

    /*!\brief Assignment via a value of a subtype that is assignable to one of the components.
     * \tparam indirect_component_type Type that models the seqan3::is_assignable_concept for
     *                                 one of the component types.
     * \param  alph                    The value of a component that should be assigned.
     *
     * ```cpp
     *     cartesian_composition<dna4, aa27> letter1{dna4::T, aa27::K};
     *
     *     letter1 = rna4::C; // yields {dna4::C, aa27::K}
     * ```
     */
    template <typename indirect_component_type>
    //!\cond
        requires one_component_is<assignable_from, indirect_component_type>
    //!\endcond
    constexpr derived_type & operator=(indirect_component_type const alph)
    {
        using component_type = meta::front<meta::find_if<components, assignable_from<indirect_component_type>>>;
        get<component_type>(*this) = alph; // delegate assignment
        return static_cast<derived_type &>(*this);
    }
    //!\cond
    // If not assignable but implicit convertible, convert first and assign afterwards
    template <typename indirect_component_type>
    //!\cond
        requires !one_component_is<assignable_from, indirect_component_type> &&
                 one_component_is<implicitely_convertible_from, indirect_component_type>
    //!\endcond
    constexpr derived_type & operator=(indirect_component_type const alph)
    {
        using component_type = meta::front<meta::find_if<components, implicitely_convertible_from<indirect_component_type>>>;
        component_type tmp(alph);
        get<component_type>(*this) = tmp;
        return static_cast<derived_type &>(*this);
    }
    //!\endcond
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the letter combination's numeric value (or "rank") in the alphabet composition.
     * \par Complexity
     * Linear in the number of alphabets.
     */
    constexpr rank_type to_rank() const
    {
        return to_rank_impl(positions);
    }

    /*!\brief Explicit cast to a single letter. Works only if the type is unique in the type list.
     * \par Complexity
     * Linear in the number of alphabets.
     */
    template <typename type>
    constexpr explicit operator type() const
    //!\cond
        requires is_unique_component<type>
    //!\endcond
    {
        return get<type>(*this);
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from a numeric value.
     * \par Complexity
     * Linear in the number of alphabets.
     * \par Exceptions
     * Asserts that the parameter is smaller than value_size [only in debug mode].
     */
    constexpr derived_type & assign_rank(rank_type const i)
    {
        assert(i < value_size);
        assign_rank_impl<0>(i);
        return static_cast<derived_type &>(*this);
    }
    //!\}

    //!\brief befriend the derived type so that it can instantiate
    //!\sa https://isocpp.org/blog/2017/04/quick-q-prevent-user-from-derive-from-incorrect-crtp-base
    friend derived_type;

    // Inherit default comparators explicitly from base class
    using base_t::operator==;
    using base_t::operator!=;
    using base_t::operator<;
    using base_t::operator>;
    using base_t::operator<=;
    using base_t::operator>=;

    /*!\name Member comparison operators with one of its component type or subtype.
     * \brief This enables the comparison of a cartesian composition value against
     * a value of the inner **alphabet type** while ignoring the other values.
     * The operators delegate to the comparison operators of the alphabet type
     * by retrieving the comparable element via `std::get<component_type>`.
     * It either takes a type that is itself a component type contained the
     * cartesian composition, or a subtype that is comparable to one of the
     * component types.
     * \{
     */
    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr bool operator==(component_type const & rhs) const noexcept
    {
        return std::get<component_type>(*this) == rhs;
    }

    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr bool operator!=(component_type const & rhs) const noexcept
    {
        return std::get<component_type>(*this) != rhs;
    }

    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr bool operator<(component_type const & rhs) const noexcept
    {
        return std::get<component_type>(*this) < rhs;
    }

    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr bool operator>(component_type const & rhs) const noexcept
    {
        return std::get<component_type>(*this) > rhs;
    }

    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr bool operator<=(component_type const & rhs) const noexcept
    {
        return std::get<component_type>(*this) <= rhs;
    }

    template <typename component_type>
    //!\cond
        requires is_unique_component<component_type>
    //!\endcond
    constexpr bool operator>=(component_type const & rhs) const noexcept
    {
        return std::get<component_type>(*this) >= rhs;
    }

    template <typename indirect_component_type>
        requires one_component_is<weakly_equality_comparable_with, indirect_component_type>
    constexpr bool operator==(indirect_component_type const rhs) const noexcept
    {
        using component_type = meta::front<meta::find_if<components, weakly_equality_comparable_with<indirect_component_type>>>;
        return get<component_type>(*this) == rhs;
    }

    template <typename indirect_component_type>
        requires one_component_is<weakly_equality_comparable_with, indirect_component_type>
    constexpr bool operator!=(indirect_component_type const & rhs) const noexcept
    {
        using component_type = meta::front<meta::find_if<components, weakly_equality_comparable_with<indirect_component_type>>>;
        return get<component_type>(*this) != rhs;
    }

    template <typename indirect_component_type>
        requires one_component_is<weakly_ordered_with, indirect_component_type>
    constexpr bool operator<(indirect_component_type const & rhs) const noexcept
    {
        using component_type = meta::front<meta::find_if<components, weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) < rhs;
    }

    template <typename indirect_component_type>
        requires one_component_is<weakly_ordered_with, indirect_component_type>
    constexpr bool operator>(indirect_component_type const & rhs) const noexcept
    {
        using component_type = meta::front<meta::find_if<components, weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) > rhs;
    }

    template <typename indirect_component_type>
        requires one_component_is<weakly_ordered_with, indirect_component_type>
    constexpr bool operator<=(indirect_component_type const & rhs) const noexcept
    {
        using component_type = meta::front<meta::find_if<components, weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) <= rhs;
    }

    template <typename indirect_component_type>
        requires one_component_is<weakly_ordered_with, indirect_component_type>
    constexpr bool operator>=(indirect_component_type const & rhs) const noexcept
    {
        using component_type = meta::front<meta::find_if<components, weakly_ordered_with<indirect_component_type>>>;
        return get<component_type>(*this) >= rhs;
    }
    //!\}

private:

    //!\brief the cumulative alphabet size products (first, first*second, first*second*third...) are cached
    static constexpr std::array<rank_type, sizeof...(component_types) + 1> cummulative_alph_sizes
    {
        [] () constexpr
        {
            std::array<rank_type, sizeof...(component_types) + 1> ret{};
            size_t count = 0;
            meta::for_each(meta::list<first_component_type, component_types...>{}, [&] (auto && alph) constexpr
            {
                ret[count] = static_cast<rank_type>(
                    alphabet_size_v<std::decay_t<decltype(alph)>> * (count > 0 ? ret[count - 1] : 1));
                ++count;
            });

            return std::move(ret);
        }()
    };

    //!\brief An index sequence up to the number of contained letters.
    static constexpr auto positions = std::make_index_sequence<sizeof...(component_types)>{};

    //!\brief Implementation of to_rank().
    template <std::size_t ...idx>
    constexpr rank_type to_rank_impl(std::index_sequence<idx...> const &) const
    {
        using seqan3::to_rank;
        if constexpr (sizeof...(idx) > 0)
        {
            return static_cast<rank_type>(
                to_rank(get<0>(*this)) +
                ((to_rank(get<idx + 1>(*this)) * cummulative_alph_sizes[idx]) + ...));
        }
        else
        {
            return to_rank(get<0>(*this));
        }
    }

    //!\brief Implementation of assign_rank().
    template <std::size_t j>
    constexpr void
    assign_rank_impl(rank_type const i)
    {
        using seqan3::assign_rank;
        if constexpr (j == 0)
        {
            assign_rank(get<j>(*this),
                          i % alphabet_size_v<meta::at_c<meta::list<first_component_type, component_types...>, j>>);
        } else
        {
            assign_rank(get<j>(*this),
                          (i / cummulative_alph_sizes[j - 1]) %
                          alphabet_size_v<meta::at_c<meta::list<first_component_type, component_types...>, j>>);
        }

        if constexpr (j < sizeof...(component_types))
            assign_rank_impl<j + 1>(i);
    }
};

/*!\name Friend comparison operators with one of its component alphabet types
 * \{
 * \brief Enables the comparison of a cartesian composition type on the right
 * hand side by delegating to the member function.
 */
template <typename component_type, typename derived_type, typename first_component_type, typename ...component_types>
    requires !std::Same<component_type, cartesian_composition<derived_type, first_component_type, component_types...>> &&
             (std::detail::WeaklyEqualityComparableWith<component_type, first_component_type> ||
             (std::detail::WeaklyEqualityComparableWith<component_type, component_types> || ...))
constexpr bool operator==(component_type const & lhs,
                          cartesian_composition<derived_type, first_component_type, component_types...> const & rhs)
{
    return rhs == lhs;
}

template <typename component_type, typename derived_type, typename first_component_type, typename ...component_types>
    requires !std::Same<component_type, cartesian_composition<derived_type, first_component_type, component_types...>> &&
             (std::detail::WeaklyEqualityComparableWith<component_type, first_component_type> ||
             (std::detail::WeaklyEqualityComparableWith<component_type, component_types> || ...))
constexpr bool operator!=(component_type const & lhs,
                          cartesian_composition<derived_type, first_component_type, component_types...> const & rhs)
{
    return rhs != lhs;
}

template <typename component_type, typename derived_type, typename first_component_type, typename ...component_types>
    requires !std::Same<component_type, cartesian_composition<derived_type, first_component_type, component_types...>> &&
             (weakly_ordered_with_concept<component_type, first_component_type> ||
             (weakly_ordered_with_concept<component_type, component_types> || ...))
constexpr bool operator<(component_type const & lhs,
                         cartesian_composition<derived_type, first_component_type, component_types...> const & rhs)
{
    return rhs > lhs;
}

template <typename component_type, typename derived_type, typename first_component_type, typename ...component_types>
    requires !std::Same<component_type, cartesian_composition<derived_type, first_component_type, component_types...>> &&
             (weakly_ordered_with_concept<component_type, first_component_type> ||
             (weakly_ordered_with_concept<component_type, component_types> || ...))
constexpr bool operator>(component_type const & lhs,
                         cartesian_composition<derived_type, first_component_type, component_types...> const & rhs)
{
    return rhs < lhs;
}

template <typename component_type, typename derived_type, typename first_component_type, typename ...component_types>
    requires !std::Same<component_type, cartesian_composition<derived_type, first_component_type, component_types...>> &&
             (weakly_ordered_with_concept<component_type, first_component_type> ||
             (weakly_ordered_with_concept<component_type, component_types> || ...))
constexpr bool operator<=(component_type const & lhs,
                          cartesian_composition<derived_type, first_component_type, component_types...> const & rhs)
{
    return rhs >= lhs;
}

template <typename component_type, typename derived_type, typename first_component_type, typename ...component_types>
    requires !std::Same<component_type, cartesian_composition<derived_type, first_component_type, component_types...>> &&
             (weakly_ordered_with_concept<component_type, first_component_type> ||
             (weakly_ordered_with_concept<component_type, component_types> || ...))
constexpr bool operator>=(component_type const & lhs,
                          cartesian_composition<derived_type, first_component_type, component_types...> const & rhs)
{
    return rhs <= lhs;
}
//!\}

} // namespace seqan3
