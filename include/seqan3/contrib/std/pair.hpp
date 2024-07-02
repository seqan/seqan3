// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::stl::pair.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#ifndef SEQAN_STD_PAIR
#define SEQAN_STD_PAIR

#include <utility>

#ifdef __cpp_lib_tuple_like

namespace seqan::stl
{

using std::pair;

} // namespace seqan::stl

#else

namespace seqan::stl
{

template <class T1, class T2>
struct pair : public std::pair<T1, T2>
{
private:
    //!\brief The underlying std::pair type.
    using base_t = std::pair<T1, T2>;

public:
    /*!\name Default constructors and assignments.
     * \{
     */
    pair() = default;                         //!< Defaulted.
    pair(pair const &) = default;             //!< Defaulted.
    pair & operator=(pair const &) = default; //!< Defaulted.
    ~pair() = default;                        //!< Defaulted.
    //!\}

    using base_t::first;
    using base_t::second;

    /*!\name Construct from arguments.
     * \{
     */
    //!@{ Constructs from arguments.
    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 &> && std::is_constructible_v<T2, U2 &>)
    constexpr pair(U1 & first, U2 & second) : base_t(first, second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr pair(U1 const & first, U2 const & second) : base_t(first, second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1> && std::is_constructible_v<T2, U2>)
    constexpr pair(U1 && first, U2 && second) : base_t(std::forward<U1>(first), std::forward<U2>(second))
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr pair(U1 const && first, U2 const && second) :
        base_t(std::forward<U1 const>(first), std::forward<U2 const>(second))
    {}
    //!@}
    //!\}

    /*!\name Construct from pair.
     * \{
     */
    //!@{ Constructs from pair.
    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 &> && std::is_constructible_v<T2, U2 &>)
    constexpr pair(pair<U1, U2> & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr pair(pair<U1, U2> const & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1> && std::is_constructible_v<T2, U2>)
    constexpr pair(pair<U1, U2> && other) : base_t(std::forward<U1>(other.first), std::forward<U2>(other.second))
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr pair(pair<U1, U2> const && other) :
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
    constexpr pair(std::pair<U1, U2> & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr pair(std::pair<U1, U2> const & other) : base_t(other.first, other.second)
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1> && std::is_constructible_v<T2, U2>)
    constexpr pair(std::pair<U1, U2> && other) : base_t(std::forward<U1>(other.first), std::forward<U2>(other.second))
    {}

    template <class U1, class U2>
        requires (std::is_constructible_v<T1, U1 const> && std::is_constructible_v<T2, U2 const>)
    constexpr pair(std::pair<U1, U2> const && other) :
        base_t(std::forward<U1 const>(other.first), std::forward<U2 const>(other.second))
    {}
    //!@}
    //!\}

    /*!\name Assign from pair.
     * \{
     */
    //!@{ Assigns from pair.
    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 &> && std::is_assignable_v<T2, U2 &>)
    constexpr pair & operator=(pair<U1, U2> & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr pair & operator=(pair<U1, U2> const & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1> && std::is_assignable_v<T2, U2>)
    constexpr pair & operator=(pair<U1, U2> && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr pair & operator=(pair<U1, U2> const && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 &> && std::is_assignable_v<T2 const, U2 &>)
    constexpr pair const & operator=(pair<U1, U2> & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr pair const & operator=(pair<U1, U2> const & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1> && std::is_assignable_v<T2 const, U2>)
    constexpr pair const & operator=(pair<U1, U2> && other) const
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr pair const & operator=(pair<U1, U2> const && other) const
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
    constexpr pair & operator=(std::pair<U1, U2> & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr pair & operator=(std::pair<U1, U2> const & other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1> && std::is_assignable_v<T2, U2>)
    constexpr pair & operator=(std::pair<U1, U2> && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1, U1 const> && std::is_assignable_v<T2, U2 const>)
    constexpr pair & operator=(std::pair<U1, U2> const && other)
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 &> && std::is_assignable_v<T2 const, U2 &>)
    constexpr pair const & operator=(std::pair<U1, U2> & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr pair const & operator=(std::pair<U1, U2> const & other) const
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1> && std::is_assignable_v<T2 const, U2>)
    constexpr pair const & operator=(std::pair<U1, U2> && other) const
    {
        first = std::move(other.first);
        second = std::move(other.second);
        return *this;
    }

    template <class U1, class U2>
        requires (std::is_assignable_v<T1 const, U1 const> && std::is_assignable_v<T2 const, U2 const>)
    constexpr pair const & operator=(std::pair<U1, U2> const && other) const
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

    /*!\name Comparison operators (pair)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator==(pair const & lhs, pair<U1, U2> const & rhs)
    {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator!=(pair const & lhs, pair<U1, U2> const & rhs)
    {
        return lhs.first != rhs.first && lhs.second != rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<(pair const & lhs, pair<U1, U2> const & rhs)
    {
        return lhs.first < rhs.first && lhs.second < rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<=(pair const & lhs, pair<U1, U2> const & rhs)
    {
        return lhs.first <= rhs.first && lhs.second <= rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>(pair const & lhs, pair<U1, U2> const & rhs)
    {
        return lhs.first > rhs.first && lhs.second > rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>=(pair const & lhs, pair<U1, U2> const & rhs)
    {
        return lhs.first >= rhs.first && lhs.second >= rhs.second;
    }

#    ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A pair with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::three_way_comparable_with<T1, U1> && std::three_way_comparable_with<T2, U2>)
    constexpr friend std::common_comparison_category_t<std::compare_three_way_result_t<U1, T1>,
                                                       std::compare_three_way_result_t<U2, T2>>
    operator<=>(pair const & lhs, pair<U1, U2> const & rhs)
    {
        if (auto cmp = lhs.first <=> rhs.first; cmp != 0)
            return cmp;
        return lhs.second <=> rhs.second;
    }
#    endif // __cpp_lib_three_way_comparison
    //!\}

    /*!\name Comparison operators (std::pair)
     * \{
     */
    /*!\brief Checks whether `lhs` and `rhs` are equal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator==(pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }

    /*!\brief Checks whether `lhs` and `rhs` are unequal.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::equality_comparable_with<T1, U1> && std::equality_comparable_with<T2, U2>)
    constexpr friend bool operator!=(pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first != rhs.first && lhs.second != rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<(pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first < rhs.first && lhs.second < rhs.second;
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator<=(pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first <= rhs.first && lhs.second <= rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>(pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first > rhs.first && lhs.second > rhs.second;
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns A bool indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::totally_ordered_with<T1, U1> && std::totally_ordered_with<T2, U2>)
    constexpr friend bool operator>=(pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        return lhs.first >= rhs.first && lhs.second >= rhs.second;
    }

#    ifdef __cpp_lib_three_way_comparison
    /*!\brief Performs a three-way comparison between `lhs` and `rhs`
     * \tparam U1 The type of the first element of `rhs`. Automatically deduced.
     * \tparam U2 The type of the second element of `rhs`. Automatically deduced.
     * \param lhs A pair.
     * \param rhs A std::pair with possibly different element types.
     * \returns An [ordering](https://en.cppreference.com/w/cpp/language/operator_comparison#Three-way_comparison)
     *          indicating the result of the comparison.
     */
    template <class U1, class U2>
        requires (std::three_way_comparable_with<T1, U1> && std::three_way_comparable_with<T2, U2>)
    constexpr friend std::common_comparison_category_t<std::compare_three_way_result_t<U1, T1>,
                                                       std::compare_three_way_result_t<U2, T2>>
    operator<=>(pair const & lhs, std::pair<U1, U2> const & rhs)
    {
        if (auto cmp = lhs.first <=> rhs.first; cmp != 0)
            return cmp;
        return lhs.second <=> rhs.second;
    }
#    endif // __cpp_lib_three_way_comparison
    //!\}
};

//!\brief Class template argument deduction guide.
template <class T1, class T2>
pair(T1, T2) -> pair<T1, T2>;

} // namespace seqan::stl

namespace std
{

template <class T1, class T2>
struct tuple_size<seqan::stl::pair<T1, T2>> : public tuple_size<std::pair<T1, T2>>
{};

template <size_t index, class T1, class T2>
struct tuple_element<index, seqan::stl::pair<T1, T2>> : public tuple_element<index, std::pair<T1, T2>>
{};

template <class T1, class T2, class U1, class U2>
    requires requires { typename seqan::stl::pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>; }
struct common_type<seqan::stl::pair<T1, T2>, seqan::stl::pair<U1, U2>>
{
    using type = seqan::stl::pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>;
};

template <class T1, class T2, class U1, class U2>
    requires requires { typename seqan::stl::pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>; }
struct common_type<std::pair<T1, T2>, seqan::stl::pair<U1, U2>>
{
    using type = seqan::stl::pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>;
};

template <class T1, class T2, class U1, class U2>
    requires requires { typename seqan::stl::pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>; }
struct common_type<seqan::stl::pair<T1, T2>, std::pair<U1, U2>>
{
    using type = seqan::stl::pair<std::common_type_t<T1, U1>, std::common_type_t<T2, U2>>;
};

template <class T1, class T2, class U1, class U2, template <class> class TQual, template <class> class UQual>
    requires requires {
        typename seqan::stl::pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                  std::common_reference_t<TQual<T2>, UQual<U2>>>;
    }
struct basic_common_reference<seqan::stl::pair<T1, T2>, seqan::stl::pair<U1, U2>, TQual, UQual>
{
    using type =
        seqan::stl::pair<std::common_reference_t<TQual<T1>, UQual<U1>>, std::common_reference_t<TQual<T2>, UQual<U2>>>;
};

template <class T1, class T2, class U1, class U2, template <class> class TQual, template <class> class UQual>
    requires requires {
        typename seqan::stl::pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                  std::common_reference_t<TQual<T2>, UQual<U2>>>;
    }
struct basic_common_reference<std::pair<T1, T2>, seqan::stl::pair<U1, U2>, TQual, UQual>
{
    using type =
        seqan::stl::pair<std::common_reference_t<TQual<T1>, UQual<U1>>, std::common_reference_t<TQual<T2>, UQual<U2>>>;
};

template <class T1, class T2, class U1, class U2, template <class> class TQual, template <class> class UQual>
    requires requires {
        typename seqan::stl::pair<std::common_reference_t<TQual<T1>, UQual<U1>>,
                                  std::common_reference_t<TQual<T2>, UQual<U2>>>;
    }
struct basic_common_reference<seqan::stl::pair<T1, T2>, std::pair<U1, U2>, TQual, UQual>
{
    using type =
        seqan::stl::pair<std::common_reference_t<TQual<T1>, UQual<U1>>, std::common_reference_t<TQual<T2>, UQual<U2>>>;
};

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan::stl::pair<T1, T2>> & get(seqan::stl::pair<T1, T2> & t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> &>(t));
}

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan::stl::pair<T1, T2>> const & get(seqan::stl::pair<T1, T2> const & t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> const &>(t));
}

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan::stl::pair<T1, T2>> && get(seqan::stl::pair<T1, T2> && t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> &&>(std::move(t)));
}

template <std::size_t i, class T1, class T2>
constexpr std::tuple_element_t<i, seqan::stl::pair<T1, T2>> const && get(seqan::stl::pair<T1, T2> const && t) noexcept
    requires (i < 2)
{
    return std::get<i>(static_cast<std::pair<T1, T2> const &&>(std::move(t)));
}

template <typename type, class T1, class T2>
constexpr type & get(seqan::stl::pair<T1, T2> & t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> &>(t));
}

template <typename type, class T1, class T2>
constexpr type const & get(seqan::stl::pair<T1, T2> const & t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> const &>(t));
}

template <typename type, class T1, class T2>
constexpr type && get(seqan::stl::pair<T1, T2> && t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> &&>(std::move(t)));
}

template <typename type, class T1, class T2>
constexpr type const && get(seqan::stl::pair<T1, T2> const && t) noexcept
    requires (!std::same_as<T1, T2>)
{
    return std::get<type>(static_cast<std::pair<T1, T2> const &&>(std::move(t)));
}

} // namespace std

#endif // ifdef __cpp_lib_tuple_like

#endif // SEQAN_STD_PAIR
