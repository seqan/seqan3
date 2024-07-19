// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Adaptations of concepts from the standard library.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <initializer_list>
#include <iterator>
#include <type_traits>

#include <seqan3/core/platform.hpp>

#if SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
#    include <string>

namespace seqan3::detail
{

//!\brief Returns whether `basic_string_t` is of type `std::basic_string<value_t, traits_t, allocator_t>`.
//!\ingroup utility_container
//!\attention Will be deleted once seqan3::detail::sequence_container_modified_by_const_iterator_bug is fixed.
template <typename basic_string_t>
struct is_basic_string : std::false_type
{};

//!\brief Returns whether `basic_string_t` is of type `std::basic_string<value_t, traits_t, allocator_t>`.
//!\ingroup utility_container
//!\attention Will be deleted once seqan3::detail::sequence_container_modified_by_const_iterator_bug is fixed.
template <typename value_t, typename traits_t, typename allocator_t>
struct is_basic_string<std::basic_string<value_t, traits_t, allocator_t>> : std::true_type
{};

//!\brief Shorthand of seqan3::detail::is_basic_string
//!\ingroup utility_container
//!\attention Will be deleted once seqan3::detail::sequence_container_modified_by_const_iterator_bug is fixed.
template <typename basic_string_t>
constexpr bool is_basic_string_v = is_basic_string<basic_string_t>::value;

} // namespace seqan3::detail
#endif // SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI

namespace seqan3
{

/*!\interface seqan3::container <>
 * \ingroup utility_container
 * \extends std::ranges::forward_range
 * \extends std::ranges::sized_range
 * \extends std::ranges::common_range
 * \extends seqan3::const_iterable_range
 * \brief The (most general) container concept as defined by the standard library.
 * \details
 * The container concept is modelled as in the [STL](https://en.cppreference.com/w/cpp/named_req/Container), but
 * additionally requires std::swap() to be implemented.
 *
 * \attention
 * Other than one might expect, `std::forward_list` does not satisfy this concept (because it does not provide
 * `.size()`).
 *
 * \noapi{Exposition only.}
 */
//!\cond
template <typename type>
concept container = requires (type val, type val2, type const cval, typename type::iterator it) {
    // member types
    typename type::value_type;
    typename type::reference;
    typename type::const_reference;
    /*
    typename type::iterator;
    requires std::forward_iterator<typename type::iterator>;
    // NOTE check whether iterator is const convertible
    {it} -> std::same_as<typename type::const_iterator>;

    typename type::const_iterator;
    requires std::forward_iterator<typename type::const_iterator>;

    typename type::difference_type;
    typename type::size_type;
    requires std::is_same_v<
        typename type::difference_type,
        typename std::iterator_traits<typename type::iterator>::difference_type
    >;
    requires std::is_same_v<
        typename std::iterator_traits<typename type::iterator>::difference_type,
        typename std::iterator_traits<typename type::const_iterator>::difference_type
    >;
*/
    // methods and operator
    { type{} } -> std::same_as<type>;       // default constructor
    { type{type{}} } -> std::same_as<type>; // copy/move constructor
    { val = val2 } -> std::same_as<type &>; // assignment
    { (&val)->~type() };                    // destructor

    { val.begin() } -> std::same_as<typename type::iterator>;
    { val.end() } -> std::same_as<typename type::iterator>;
    { cval.begin() } -> std::same_as<typename type::const_iterator>;
    { cval.end() } -> std::same_as<typename type::const_iterator>;
    { val.cbegin() } -> std::same_as<typename type::const_iterator>;
    { val.cend() } -> std::same_as<typename type::const_iterator>;

    requires !std::equality_comparable<typename type::value_type> || std::equality_comparable<type>;

    { val.swap(val2) } -> std::same_as<void>;
    { swap(val, val2) } -> std::same_as<void>;
    { std::swap(val, val2) } -> std::same_as<void>;

    { val.size() } -> std::same_as<typename type::size_type>;
    { val.max_size() } -> std::same_as<typename type::size_type>;
    { val.empty() } -> std::same_as<bool>;
};
//!\endcond

/*!\interface seqan3::sequence_container <>
 * \ingroup utility_container
 * \extends seqan3::container
 * \brief A more refined container concept than seqan3::container.
 *
 * Includes constraints on constructors, `assign()`, `.insert()`, `.erase()`, `.push_back()`, `.pop_back`, `.clear()`,
 * `.size()`, `front()` and `.back()` member functions with corresponding signatures. Models the subset of the
 * [STL SequenceConcept](https://en.cppreference.com/w/cpp/named_req/SequenceContainer) that is supported
 * by `std::list`, `std::vector`, `std::deque`, `std::basic_string`.
 *
 * \attention
 * `std::array` and `std::forward_list` do not satisfy this concept.
 *
 * \noapi{Exposition only.}
 */
//!\cond
template <typename type>
concept sequence_container = requires (type val, type val2, type const cval) {
    requires container<type>;

    // construction
    { type(typename type::size_type{}, typename type::value_type{}) };
    { type{val2.begin(), val2.end()} }; // NOTE that this could be any input iterator:
    { type{std::initializer_list<typename type::value_type>{}} };
    { val = std::initializer_list<typename type::value_type>{} } -> std::same_as<type &>;

    // assignment NOTE return type is type & for std::string and void for other stl containers:
    { val.assign(val2.begin(), val2.end()) };
    { val.assign(std::initializer_list<typename type::value_type>{}) };
    { val.assign(typename type::size_type{}, typename type::value_type{}) };

    // modify container
    // TODO: how do you model this?
    // {val.emplace(typename type::const_iterator{}, ?)} -> std::same_as<typename type::iterator>;
#if SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
    { val.insert(val.begin(), val2.front()) } -> std::same_as<typename type::iterator>;
    { val.insert(val.begin(), typename type::value_type{}) } -> std::same_as<typename type::iterator>;

    // std::string is missing the const_iterator versions for insert in pre-C++11 ABI
    requires detail::is_basic_string_v<type>
                 || requires (type val, type val2)
#else  // ^^^ workaround / no workaround vvv
    requires requires (type val, type val2)
#endif // SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
    {
        { val.insert(val.cbegin(), val2.front()) } -> std::same_as<typename type::iterator>;
        { val.insert(val.cbegin(), typename type::value_type{}) } -> std::same_as<typename type::iterator>;
        {
            val.insert(val.cbegin(), typename type::size_type{}, typename type::value_type{})
        } -> std::same_as<typename type::iterator>;
        { val.insert(val.cbegin(), val2.begin(), val2.end()) } -> std::same_as<typename type::iterator>;
        {
            val.insert(val.cbegin(), std::initializer_list<typename type::value_type>{})
        } -> std::same_as<typename type::iterator>;
    };

#if SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
    // std::string is missing the const_iterator versions for erase in pre-C++11 ABI
    requires detail::is_basic_string_v<type>
                 || requires (type val)
#else  // ^^^ workaround / no workaround vvv
    requires requires (type val)
#endif // SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
    {
        { val.erase(val.cbegin()) } -> std::same_as<typename type::iterator>;
        { val.erase(val.cbegin(), val.cend()) } -> std::same_as<typename type::iterator>;
    };

    { val.push_back(val.front()) } -> std::same_as<void>;
    { val.push_back(typename type::value_type{}) } -> std::same_as<void>;
    { val.pop_back() } -> std::same_as<void>;
    { val.clear() } -> std::same_as<void>;

    // access container
    { val.front() } -> std::same_as<typename type::reference>;
    { cval.front() } -> std::same_as<typename type::const_reference>;
    { val.back() } -> std::same_as<typename type::reference>;
    { cval.back() } -> std::same_as<typename type::const_reference>;
};
//!\endcond

/*!\interface seqan3::random_access_container <>
 * \ingroup utility_container
 * \extends seqan3::sequence_container
 * \extends std::ranges::random_access_range
 * \brief A more refined container concept than seqan3::sequence_container.
 *
 * Adds requirements for `.at()`, `.resize()` and the subscript operator `[]`. Models the subset of the
 * [STL SequenceConcept](https://en.cppreference.com/w/cpp/named_req/SequenceContainer) that is supported
 * by `std::vector`, `std::deque` and `std::basic_string`.
 *
 * \attention
 * `std::array`, `std::forward_list` and `std::list` do not satisfy this concept.
 *
 * \noapi{Exposition only.}
 */
//!\cond
template <typename type>
concept random_access_container = requires (type val) {
    requires sequence_container<type>;

    // access container
    { val[0] } -> std::same_as<typename type::reference>;
    { val.at(0) } -> std::same_as<typename type::reference>;

    // modify container
    { val.resize(0) } -> std::same_as<void>;
    { val.resize(0, typename type::value_type{}) } -> std::same_as<void>;
};
//!\endcond

/*!\interface seqan3::reservible_container <>
 * \ingroup utility_container
 * \extends seqan3::random_access_container
 * \brief A more refined container concept than seqan3::random_access_container.
 *
 * Adds requirements for `.reserve()`, `.capacity()` and `.shrink_to_fit()`.
 * Satisfied by `std::vector` and `std::basic_string`.
 *
 * \attention
 * `std::array`, `std::forward_list`, `std::list` and `std::deque` do not satisfy this concept.
 *
 * \noapi{Exposition only.}
 */
//!\cond
template <typename type>
concept reservible_container = requires (type val) {
    requires random_access_container<type>;

    { val.capacity() } -> std::same_as<typename type::size_type>;
    { val.reserve(0) } -> std::same_as<void>;
    { val.shrink_to_fit() } -> std::same_as<void>;
};
//!\endcond

} // namespace seqan3
