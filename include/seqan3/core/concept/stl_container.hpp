// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

#pragma once

#include <initializer_list>

/*! \file core/concept/stl_container.hpp
 * \brief Adaptations of concepts from the standard library
 * \ingroup core
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

namespace seqan3
{

/*! \var concept bool container_concept
 *  \brief The (most general) container concept
 *
 * The container concept is modelled exactly as in the [STL](http://en.cppreference.com/w/cpp/concept/Container).
 * \privatesection
 */

template <typename type>
concept bool container_concept = requires (type val, type val2)
{
    // member types
    typename type::value_type;
    typename type::reference;
    typename type::const_reference;
    typename type::iterator; //TODO must satisfy forward_iterator_concept and convertible to const_interator
    typename type::const_iterator; //TODO must satisfy forward_iterator_concept
    typename type::difference_type;
    typename type::size_type; // TODO must be the same as iterator_traits::difference_type for iterator and const_iterator

    // methods and operator
    { type{}          } -> type;   // default constructor
    { type{type{}}    } -> type;   // copy/move constructor
    { val = val2      } -> type &; // assignment
    { (&val)->~type() } -> void;   // destructor

    { val.begin()     } -> typename type::iterator;
    { val.end()       } -> typename type::iterator;
    { val.cbegin()    } -> typename type::const_iterator;
    { val.cend()      } -> typename type::const_iterator;

    { val == val2     } -> bool;
    { val != val2     } -> bool;

    { val.swap(val2)  } -> void;
    { swap(val, val2) } -> void;

    { val.max_size()  } -> typename type::size_type;
    { val.empty()     } -> bool;
};

/*!
 * \publicsection
 * \var concept bool sequence_light_concept
 * \brief A more refined container concept than @link container_concept @endlink
 *
 * Includes contraints on constructors, assign member and `.front()`. Models the subset of the
 * [STL SequenceContainerConcept](http://en.cppreference.com/w/cpp/concept/SequenceContainer) that is supported
 * by `std::forward_list`, `std::list`, `std::vector`, `std::deque`, `std::basic_string`.
 *
 * `std::array` does not satisfy this concept.
 * \privatesection
 */

template <typename type>
concept bool sequence_light_concept = requires (type val, type val2)
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

    { val.front() } -> typename type::value_type &;
};

/*!
 * \publicsection
 * \var concept bool sequence_concept
 * \brief A more refined container concept than @link sequence_light_concept @endlink
 *
 * Adds requirements for `.insert()`, `.erase()`, `.push_back()`, `.pop_back`, `.clear()`, `.size()` and `.back()`
 * member functions with corresponding signatures. Models the subset of the
 * [STL SequenceContainerConcept](http://en.cppreference.com/w/cpp/concept/SequenceContainer) that is supported
 * by `std::list`, `std::vector`, `std::deque` and `std::basic_string`.
 *
 * `std::array` and `std::forward_list` do not satisfy this concept.
 * \privatesection
 */

template <typename type>
concept bool sequence_concept = requires (type val, type val2)
{
    requires sequence_light_concept<type>;

    // modify container
//TODO: how do you model this?
//     { val.emplace(typename type::const_iterator{}, ?                                   } -> typename type::iterator;
    { val.insert(val.begin(), val2.front())                                            } -> typename type::iterator;
    { val.insert(val.begin(), typename type::value_type{})                             } -> typename type::iterator;
    { val.insert(val.cbegin(), typename type::size_type{}, typename type::value_type{})} -> typename type::iterator;
    { val.insert(val.cbegin(), val2.begin(), val2.end())                               } -> typename type::iterator;
//TODO this fails on std::string, although it should work
//     { val.insert(val.cbegin(), std::initializer_list<typename type::value_type>{})    } -> typename type::iterator;
    { val.erase(val.cbegin())                                                          } -> typename type::iterator;
    { val.erase(val.cbegin(), val.cend())                                              } -> typename type::iterator;
    { val.push_back(val.front())                                                       } -> void;
    { val.push_back(typename type::value_type{})                                       } -> void;
    { val.pop_back()                                                                   } -> void;
    { val.clear()                                                                      } -> void;

    // access container
    { val.size() } -> typename type::size_type;
    { val.back() } -> typename type::value_type &;
};

/*!
 * \publicsection
 * \var concept bool random_access_sequence_concept
 * \brief A more refined container concept than @link sequence_concept @endlink
 *
 * Adds requirements for `.at()`, `.resize()` and the subscript operator `[]`. Models the subset of the
 * [STL SequenceContainerConcept](http://en.cppreference.com/w/cpp/concept/SequenceContainer) that is supported
 * by `std::vector`, `std::deque` and `std::basic_string`.
 *
 * `std::array`, `std::forward_list` and `std::list` do not satisfy this concept.
 * \privatesection
 */

template <typename type>
concept bool random_access_sequence_concept = requires (type val)
{
    requires sequence_concept<type>;

    // access container
    { val[0]    } -> typename type::value_type &;
    { val.at(0) } -> typename type::value_type &;

    // modify container
    { val.resize(0)                              } -> void;
    { val.resize(0, typename type::value_type{}) } -> void;
};

/*!
 * \publicsection
 * \var concept bool container_of_container_concept
 * \brief A multi-dimensional @link container_concept @endlink.
 *
 * Requires that both the type and it's `value_type` fulfill @link container_concept @endlink.
 * \privatesection
 */

template <typename type>
concept bool container_of_container_concept = requires (type val)
{
    requires container_concept<type>;
    requires container_concept<typename type::value_type>;
};

/*!
 * \publicsection
 * \var concept bool sequence_of_sequence_concept
 * \brief A multi-dimensional @link sequence_concept @endlink.
 *
 * Requires that both the type and it's `value_type` fulfill @link sequence_concept @endlink.
 * \privatesection
 */

template <typename type>
concept bool sequence_of_sequence_concept = requires (type val)
{
    requires sequence_concept<type>;
    requires sequence_concept<typename type::value_type>;
};

/*!
 * \publicsection
 * \var concept bool ra_sequence_of_ra_sequence_concept
 * \brief A multi-dimensional @link random_access_sequence_concept @endlink
 *
 * Requires that both the type and it's `value_type` fulfill @link random_access_sequence_concept @endlink.
 * \privatesection
 */

template <typename type>
concept bool ra_sequence_of_ra_sequence_concept = requires (type val)
{
    requires random_access_sequence_concept<type>;
    requires random_access_sequence_concept<typename type::value_type>;
};

} // namespace seqan3

#ifndef NDEBUG
/* Check the STL containers */

#include <vector>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <string>

static_assert(seqan3::container_concept<std::array<char, 2>>);
static_assert(seqan3::sequence_light_concept<std::forward_list<char>>);
static_assert(seqan3::sequence_concept<std::list<char>>);
static_assert(seqan3::random_access_sequence_concept<std::vector<char>>);
static_assert(seqan3::random_access_sequence_concept<std::deque<char>>);
static_assert(seqan3::random_access_sequence_concept<std::string>);

static_assert(seqan3::container_of_container_concept<std::array<std::array<char, 2>, 2>>);
static_assert(seqan3::sequence_of_sequence_concept<std::list<std::list<char>>>);
static_assert(seqan3::ra_sequence_of_ra_sequence_concept<std::vector<std::vector<char>>>);

/* Check the SDSL containers */
//TODO

/* Check our containers */
//TODO
#endif
