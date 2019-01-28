// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Adaptations of concepts from the standard library.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <initializer_list>
#include <iterator>
#include <type_traits>

// remove if sequence_container_concept_modified_by_const_iterator_bug vanished from travis
#include <string>

#include <seqan3/std/iterator>

// TODO:
// * merge sequence_container_concept_modified_by_const_iterator back into
//   sequence_container_concept
// * remove is_basic_string
// * fix test cases
// * remove #include <string> in this file
// once the ubuntu::ppa [1] of g++-7 has a newer update than
// 7.2.0-1ubuntu1~16.04 (2017-08-20)
//
// [1] https://launchpad.net/~ubuntu-toolchain-r/+archive/ubuntu/test?field.series_filter=xenial
namespace seqan3::detail
{
//!\privatesection

//!\brief Returns whether `basic_string_t` is of type `std::basic_string<value_t, traits_t, allocator_t>`.
//!\attention Will be deleted once seqan3::detail::sequence_container_concept_modified_by_const_iterator_bug is fixed.
template <typename basic_string_t>
struct is_basic_string : std::false_type
{};

//!\brief Returns whether `basic_string_t` is of type `std::basic_string<value_t, traits_t, allocator_t>`.
//!\attention Will be deleted once seqan3::detail::sequence_container_concept_modified_by_const_iterator_bug is fixed.
template <typename value_t, typename traits_t, typename allocator_t>
struct is_basic_string<std::basic_string<value_t, traits_t, allocator_t>> : std::true_type
{};

//!\brief Shorthand of seqan3::detail::is_basic_string
//!\attention Will be deleted once seqan3::detail::sequence_container_concept_modified_by_const_iterator_bug is fixed.
template <typename basic_string_t>
constexpr bool is_basic_string_v = is_basic_string<basic_string_t>::value;

/*!\interface seqan3::detail::sequence_container_concept_modified_by_const_iterator <>
 * \brief Checks whether insert and erase can be used with const_iterator
 *
 * \attention This will be merged back into sequence_container_concept once
 * seqan3::detail::sequence_container_concept_modified_by_const_iterator_bug is fixed.
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT sequence_container_concept_modified_by_const_iterator = requires (type val, type val2)
{
    { val.insert(val.cbegin(), val2.front())                                           } -> typename type::iterator;
    { val.insert(val.cbegin(), typename type::value_type{})                            } -> typename type::iterator;
    { val.insert(val.cbegin(), typename type::size_type{}, typename type::value_type{})} -> typename type::iterator;
    { val.insert(val.cbegin(), val2.begin(), val2.end())                               } -> typename type::iterator;
    requires is_basic_string_v<type> || requires(type val)
    {
        // TODO this function is not defined on strings (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=83328)
        { val.insert(val.cbegin(), std::initializer_list<typename type::value_type>{}) } -> typename type::iterator;
    };
    { val.erase(val.cbegin())                                                          } -> typename type::iterator;
    { val.erase(val.cbegin(), val.cend())                                              } -> typename type::iterator;

    { val.insert(val.begin(), typename type::size_type{}, typename type::value_type{}) } -> typename type::iterator;
    { val.insert(val.begin(), val2.begin(), val2.end())                                } -> typename type::iterator;
};
//!\endcond

/*!\brief Workaround for a ubuntu/travis-ci exclusive bug with g++-7.2.
 *
 * seqan3::detail::sequence_container_concept_modified_by_const_iterator <std::string> is
 * known to work, but ubuntu::ppa (<18.04)/travis-ci has a version of g++-7.2
 * where a bug in the STL prevents this concept to be true.
 *
 * \attention This workaround can be removed if
 * `/test/range/container/container_concept_test.cpp` is not failing on
 * ubuntu::ppa (<18.04)/travis-ci anymore. \n
 * Probably when the ppa version of gcc7 is newer than `7.2.0-1ubuntu1~16.04` (2017-08-20)
 * \sa https://launchpad.net/~ubuntu-toolchain-r/+archive/ubuntu/test?field.series_filter=xenial
 */
template<typename string_t = std::string>
constexpr bool sequence_container_concept_modified_by_const_iterator_bug =
    is_basic_string_v<string_t> && !sequence_container_concept_modified_by_const_iterator<string_t>;

//!\publicsection

} // seqan3::detail

namespace seqan3
{

/*!\addtogroup container
 * \{
 */
/*!\interface seqan3::container_concept <>
 * \extends std::ranges::ForwardRange
 * \extends std::ranges::SizedRange
 * \extends std::ranges::CommonRange
 * \extends seqan3::const_iterable_concept
 * \brief The (most general) container concept as defined by the standard library.
 * \details
 * The container concept is modelled as in the [STL](http://en.cppreference.com/w/cpp/concept/Container), but
 * additionally requires std::swap() to be implemented.
 *
 * \attention
 * Other than one might expect, `std::forward_list` does not satisfy this concept (because it does not provide
 * `.size()`).
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT container_concept = requires (type val, type val2, type const cval, typename type::iterator it)
{
    // member types
    typename type::value_type;
    typename type::reference;
    typename type::const_reference;
/*
    typename type::iterator;
    requires std::ForwardIterator<typename type::iterator>;
    { it } -> typename type::const_iterator; // NOTE check whether iterator is const convertible

    typename type::const_iterator;
    requires std::ForwardIterator<typename type::const_iterator>;

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
    { type{}          } -> type;   // default constructor
    { type{type{}}    } -> type;   // copy/move constructor
    { val = val2      } -> type &; // assignment
    { (&val)->~type() } -> void;   // destructor

    { val.begin()     } -> typename type::iterator;
    { val.end()       } -> typename type::iterator;
    { cval.begin()    } -> typename type::const_iterator;
    { cval.end()      } -> typename type::const_iterator;
    { val.cbegin()    } -> typename type::const_iterator;
    { val.cend()      } -> typename type::const_iterator;

    { val == val2     } -> bool;
    { val != val2     } -> bool;

    { val.swap(val2)  } -> void;
    { swap(val, val2) } -> void;
    { std::swap(val, val2) } -> void;

    { val.size()      } -> typename type::size_type;
    { val.max_size()  } -> typename type::size_type;
    { val.empty()     } -> bool;
};
//!\endcond

/*!\interface seqan3::sequence_container_concept <>
 * \extends seqan3::container_concept
 * \brief A more refined container concept than seqan3::container_concept.
 *
 * Includes constraints on constructors, `assign()`, `.insert()`, `.erase()`, `.push_back()`, `.pop_back`, `.clear()`,
 * `.size()`, `front()` and `.back()` member functions with corresponding signatures. Models the subset of the
 * [STL SequenceContainerConcept](http://en.cppreference.com/w/cpp/concept/SequenceContainer) that is supported
 * by `std::list`, `std::vector`, `std::deque`, `std::basic_string`.
 *
 * \attention
 * `std::array` and `std::forward_list` do not satisfy this concept.
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT sequence_container_concept = requires (type val, type val2, type const cval)
{
    requires container_concept<type>;

    // construction
    { type{typename type::size_type{}, typename type::value_type{}} };
    { type{val2.begin(), val2.end()}                                }; // NOTE that this could be any input iterator:
    { type{std::initializer_list<typename type::value_type>{}}      };
    { val = std::initializer_list<typename type::value_type>{}      } -> type &;

    // assignment NOTE return type is type & for std::string and void for other stl containers:
    { val.assign(val2.begin(), val2.end())                                };
    { val.assign(std::initializer_list<typename type::value_type>{})      };
    { val.assign(typename type::size_type{}, typename type::value_type{}) };

    // modify container
    // TODO: how do you model this?
    // { val.emplace(typename type::const_iterator{}, ?                                   } -> typename type::iterator;

    { val.insert(val.begin(), val2.front())                                            } -> typename type::iterator;
    { val.insert(val.begin(), typename type::value_type{})                             } -> typename type::iterator;
    // because of a travis bug we can't assume typename type::iterator as return type
    { val.insert(val.begin(), typename type::size_type{}, typename type::value_type{}) };
    { val.insert(val.begin(), val2.begin(), val2.end())                                };
    //TODO should return type::iterator on strings (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=83328)
    { val.insert(val.begin(), std::initializer_list<typename type::value_type>{})      };
    { val.erase(val.begin())                                                           } -> typename type::iterator;
    { val.erase(val.begin(), val.end())                                                } -> typename type::iterator;

    // workaround a travis bug where insert/erase can't take a const iterator, e.g. cbegin()
    requires detail::sequence_container_concept_modified_by_const_iterator_bug<type> ||
             detail::sequence_container_concept_modified_by_const_iterator<type>;

    { val.push_back(val.front())                                                       } -> void;
    { val.push_back(typename type::value_type{})                                       } -> void;
    { val.pop_back()                                                                   } -> void;
    { val.clear()                                                                      } -> void;

    // access container
    { val.front()  } -> typename type::reference;
    { val.front()  } -> typename type::const_reference;
    { cval.front() } -> typename type::const_reference;
    { val.back()   } -> typename type::reference;
    { val.back()   } -> typename type::const_reference;
    { cval.back()  } -> typename type::const_reference;
};
//!\endcond

/*!\interface seqan3::random_access_container_concept <>
 * \extends seqan3::sequence_container_concept
 * \extends std::ranges::RandomAccessRange
 * \brief A more refined container concept than seqan3::sequence_container_concept.
 *
 * Adds requirements for `.at()`, `.resize()` and the subscript operator `[]`. Models the subset of the
 * [STL SequenceContainerConcept](http://en.cppreference.com/w/cpp/concept/SequenceContainer) that is supported
 * by `std::vector`, `std::deque` and `std::basic_string`.
 *
 * \attention
 * `std::array`, `std::forward_list` and `std::list` do not satisfy this concept.
 *
 * \sa
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT random_access_container_concept = requires (type val)
{
    requires sequence_container_concept<type>;

    // access container
    { val[0]    } -> typename type::reference;
    { val.at(0) } -> typename type::reference;

    // modify container
    { val.resize(0)                              } -> void;
    { val.resize(0, typename type::value_type{}) } -> void;
};
//!\endcond

/*!\interface seqan3::reservable_container_concept <>
 * \extends seqan3::random_access_container_concept
 * \brief A more refined container concept than seqan3::random_access_container_concept.
 *
 * Adds requirements for `.reserve()`, `.capacity()` and `.shrink_to_fit()`.
 * Satisfied by `std::vector` and `std::basic_string`.
 *
 * \attention
 * `std::array`, `std::forward_list`, `std::list` and `std::deque` do not satisfy this concept.
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT reservable_container_concept = requires (type val)
{
    requires random_access_container_concept<type>;

    { val.capacity()      } -> typename type::size_type;
    { val.reserve(0)      } -> void;
    { val.shrink_to_fit() } -> void;
};
//!\endcond

//!\}

} // namespace seqan3
