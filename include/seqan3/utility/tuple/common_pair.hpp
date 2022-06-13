// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::common_pair.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief A [std::pair](https://en.cppreference.com/w/cpp/utility/pair) implementation that incorporates most changes
 *        from C++23's standard library.
 * \ingroup utility_tuple
 */
template <class T1, class T2>
struct common_pair : public std::pair<T1, T2>
{
private:
    //!\brief The underlying std::pair type.
    using base_t = std::pair<T1, T2>;

public:
    /*!\name Default constructors and assignments.
     * \{
     */
    common_pair() = default;                                //!< Defaulted.
    common_pair(common_pair const &) = default;             //!< Defaulted.
    common_pair & operator=(common_pair const &) = default; //!< Defaulted.
    ~common_pair() = default;                               //!< Defaulted.
    //!\}

    using base_t::first;
    using base_t::second;

    /*!\name Construct from arguments.
     * \{
     */
    //!@{ Constructs from arguments.
    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 &> && std::is_constructible_v<T2, U2 &>)
    constexpr common_pair(U1 & first, U2 & second) : base_t(first, second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr common_pair(U1 const & first, U2 const & second) : base_t(first, second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1> && std::is_constructible_v<T2, U2>)
    constexpr common_pair(U1 && first, U2 && second) : base_t(std::forward<U1>(first), std::forward<U2>(second))
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr common_pair(U1 const && first, U2 const && second) :
        base_t(std::forward<U1 const>(first), std::forward<U2 const>(second))
    {}
    //!@}
    //!\}

    /*!\name Construct from common_pair.
     * \{
     */
    //!@{ Constructs from common_pair.
    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 &> && std::is_constructible_v<T2, U2 &>)
    constexpr common_pair(common_pair<U1, U2> & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr common_pair(common_pair<U1, U2> const & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1> && std::is_constructible_v<T2, U2>)
    constexpr common_pair(common_pair<U1, U2> && other) :
        base_t(std::forward<U1>(other.first), std::forward<U2>(other.second))
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr common_pair(common_pair<U1, U2> const && other) :
        base_t(std::forward<U1 const>(other.first), std::forward<U2 const>(other.second))
    {}
    //!@}
    //!\}

    /*!\name Construct from std::pair.
     * \{
     */
    //!@{ Constructs from std::pair.
    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 &> && std::is_constructible_v<T2, U2 &>)
    constexpr common_pair(std::pair<U1, U2> & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr common_pair(std::pair<U1, U2> const & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1> && std::is_constructible_v<T2, U2>)
    constexpr common_pair(std::pair<U1, U2> && other) :
        base_t(std::forward<U1>(other.first), std::forward<U2>(other.second))
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr common_pair(std::pair<U1, U2> const && other) :
        base_t(std::forward<U1 const>(other.first), std::forward<U2 const>(other.second))
    {}
    //!@}
    //!\}

    /*!\name Assign from common_pair.
     * \{
     */
    //!@{ Assigns from common_pair.
    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 &> && std::is_assignable_v<T2, U2 &>)
    constexpr common_pair & operator=(common_pair<U1, U2> & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr common_pair & operator=(common_pair<U1, U2> const & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1> && std::is_assignable_v<T2, U2>)
    constexpr common_pair & operator=(common_pair<U1, U2> && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr common_pair & operator=(common_pair<U1, U2> const && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 &> && std::is_assignable_v<T2 const, U2 &>)
    constexpr common_pair const & operator=(common_pair<U1, U2> & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr common_pair const & operator=(common_pair<U1, U2> const & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1> && std::is_assignable_v<T2 const, U2>)
    constexpr common_pair const & operator=(common_pair<U1, U2> && other) const
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr common_pair const & operator=(common_pair<U1, U2> const && other) const
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }
    //!@}
    //!\}

    /*!\name Assign from std::pair.
     * \{
     */
    //!@{ Assigns from std::pair.
    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 &> && std::is_assignable_v<T2, U2 &>)
    constexpr common_pair & operator=(std::pair<U1, U2> & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr common_pair & operator=(std::pair<U1, U2> const & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1> && std::is_assignable_v<T2, U2>)
    constexpr common_pair & operator=(std::pair<U1, U2> && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr common_pair & operator=(std::pair<U1, U2> const && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 &> && std::is_assignable_v<T2 const, U2 &>)
    constexpr common_pair const & operator=(std::pair<U1, U2> & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr common_pair const & operator=(std::pair<U1, U2> const & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1> && std::is_assignable_v<T2 const, U2>)
    constexpr common_pair const & operator=(std::pair<U1, U2> && other) const
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr common_pair const & operator=(std::pair<U1, U2> const && other) const
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }
    //!@}
    //!\}

    /*!\name Conversion to std::pair.
     * \{
     */
    //!@{ Converts to std::pair
    template <class U1, class U2>
        requires (std::is_constructible_v<U1, T1 &> && std::is_constructible_v<U2, T2 &>)
    operator std::pair<U1, U2>() &
    {
        return {first, second};
    }

    template <class U1, class U2>
        requires (std::is_constructible_v<U1, T1 const> && std::is_constructible_v<U2, T2 const>)
    operator std::pair<U1, U2>() const &
    {
        return {first, second};
    }

    template <class U1, class U2>
        requires (std::is_constructible_v<U1, T1> && std::is_constructible_v<U2, T2>)
    operator std::pair<U1, U2>() &&
    {
        return {std::move(first), std::move(second)};
    }

    template <class U1, class U2>
        requires (std::is_constructible_v<U1, T1 const> && std::is_constructible_v<U2, T2 const>)
    operator std::pair<U1, U2>() const &&
    {
        return {std::move(first), std::move(second)};
    }
    //!@}
    //!\}

    /*!\name Comparison operators (common_pair)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A common_pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator==(common_pair const & lhs, common_pair<U1, U2> const & rhs)
    {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A common_pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator!=(common_pair const & lhs, common_pair<U1, U2> const & rhs)
    {
        return lhs.first != rhs.first && lhs.second != rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A common_pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<(common_pair const & lhs, common_pair<U1, U2> const & rhs)
    {
        return lhs.first < rhs.first && lhs.second < rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A common_pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<=(common_pair const & lhs, common_pair<U1, U2> const & rhs)
    {
        return lhs.first <= rhs.first && lhs.second <= rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A common_pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>(common_pair const & lhs, common_pair<U1, U2> const & rhs)
    {
        return lhs.first > rhs.first && lhs.second > rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A common_pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>=(common_pair const & lhs, common_pair<U1, U2> const & rhs)
    {
        return lhs.first >= rhs.first && lhs.second >= rhs.second;
    }

#ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A common_pair with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::three_way_comparable_with<T1, U1> && std::three_way_comparable_with<T2, U2>)
    constexpr friend std::common_comparison_category_t<std::compare_three_way_result_t<U1, T1>,
                                                       std::compare_three_way_result_t<U2, T2>>
    operator<=>(common_pair const & lhs, common_pair<U1, U2> const & rhs)
    {
        if (auto cmp = lhs.first <=> rhs.first; cmp != 0)
            return cmp;
        return lhs.second <=> rhs.second;
    }
#endif // __cpp_lib_three_way_comparison
    //!\}

    /*!\name Comparison operators (std::pair)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator==(common_pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator!=(common_pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first != rhs.first && lhs.second != rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<(common_pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first < rhs.first && lhs.second < rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<=(common_pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first <= rhs.first && lhs.second <= rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>(common_pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first > rhs.first && lhs.second > rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>=(common_pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first >= rhs.first && lhs.second >= rhs.second;
    }

#ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A common_pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::three_way_comparable_with<T1, U1> && std::three_way_comparable_with<T2, U2>)
    constexpr friend std::common_comparison_category_t<std::compare_three_way_result_t<U1, T1>,
                                                       std::compare_three_way_result_t<U2, T2>>
    operator<=>(common_pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        if (auto cmp = lhs.first <=> rhs.first; cmp != 0)
            return cmp;
        return lhs.second <=> rhs.second;
    }
#endif // __cpp_lib_three_way_comparison
    //!\}
};

//!\brief Class template argument deduction guide.
template <class T1, class T2>
common_pair(T1, T2) -> common_pair<T1, T2>;

} // namespace seqan3

//!\cond
namespace std
{

template <class T1, class T2>
struct tuple_size<seqan3::common_pair<T1, T2>> : public tuple_size<std::pair<T1, T2>>
{};

template <size_t index, class T1, class T2>
struct tuple_element<index, seqan3::common_pair<T1, T2>> : public tuple_element<index, std::pair<T1, T2>>
{};

template <class T1, class T2, class U1, class U2>
    requires requires { typename seqan3::common_pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>; }
struct common_type<seqan3::common_pair<T1, T2>, seqan3::common_pair<U1, U2>>
{
    using type = seqan3::common_pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>;
};

template <class T1, class T2, class U1, class U2>
    requires requires { typename seqan3::common_pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>; }
struct common_type<std::pair<T1, T2>, seqan3::common_pair<U1, U2>>
{
    using type = seqan3::common_pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>;
};

template <class T1, class T2, class U1, class U2>
    requires requires { typename seqan3::common_pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>; }
struct common_type<seqan3::common_pair<T1, T2>, std::pair<U1, U2>>
{
    using type = seqan3::common_pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>;
};

template <class T1, class T2, class U1, class U2, template <class> class TQual, template <class> class UQual>
    requires requires {
                 typename seqan3::common_pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                              std::common_reference_t<TQual<T2>, UQual<U2>>>;
             }
struct basic_common_reference<seqan3::common_pair<T1, T2>, seqan3::common_pair<U1, U2>, TQual, UQual>
{
    using type = seqan3::common_pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                     std::common_reference_t<TQual<T2>, UQual<U2>>>;
};

template <class T1, class T2, class U1, class U2, template <class> class TQual, template <class> class UQual>
    requires requires {
                 typename seqan3::common_pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                              std::common_reference_t<TQual<T2>, UQual<U2>>>;
             }
struct basic_common_reference<std::pair<T1, T2>, seqan3::common_pair<U1, U2>, TQual, UQual>
{
    using type = seqan3::common_pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                     std::common_reference_t<TQual<T2>, UQual<U2>>>;
};

template <class T1, class T2, class U1, class U2, template <class> class TQual, template <class> class UQual>
    requires requires {
                 typename seqan3::common_pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                              std::common_reference_t<TQual<T2>, UQual<U2>>>;
             }
struct basic_common_reference<seqan3::common_pair<T1, T2>, std::pair<U1, U2>, TQual, UQual>
{
    using type = seqan3::common_pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                     std::common_reference_t<TQual<T2>, UQual<U2>>>;
};

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan3::common_pair<T1, T2>> & get(seqan3::common_pair<T1, T2> & t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> &>(t));
}

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan3::common_pair<T1, T2>> const &
get(seqan3::common_pair<T1, T2> const & t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> const &>(t));
}

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan3::common_pair<T1, T2>> && get(seqan3::common_pair<T1, T2> && t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> &&>(std::move(t)));
}

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan3::common_pair<T1, T2>> const &&
get(seqan3::common_pair<T1, T2> const && t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> const &&>(std::move(t)));
}

template <typename type, class T1, class T2>
constexpr type & get(seqan3::common_pair<T1, T2> & t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> &>(t));
}

template <typename type, class T1, class T2>
constexpr type const & get(seqan3::common_pair<T1, T2> const & t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> const &>(t));
}

template <typename type, class T1, class T2>
constexpr type && get(seqan3::common_pair<T1, T2> && t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> &&>(std::move(t)));
}

template <typename type, class T1, class T2>
constexpr type const && get(seqan3::common_pair<T1, T2> const && t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> const &&>(std::move(t)));
}

} // namespace std
//!\endcond
