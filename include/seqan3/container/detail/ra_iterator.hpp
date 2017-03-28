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
// Authors: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
//          Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ============================================================================

#pragma once

#include <seqan3/container/concepts.hpp>

namespace seqan3::detail
{

/*
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
*/

template <typename container_type>
    requires container_concept<container_type> //requires random_access_range_concept<container_type>
struct ra_iterator
{

private:
    container_type & host;
    typename container_type::size_type pos;
public:
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::value_type value_type;
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

    //! default constructor
    ra_iterator() : host(*(container_type*)0) {};
    //! initialize with pointer of container structure
    ra_iterator(container_type & _host) : host{_host}, pos{0} {};
    //! copy constructors
    ra_iterator(ra_iterator const &) = default;
    ra_iterator(ra_iterator &&) = default;
    ra_iterator& operator=(ra_iterator const &) = default;
    ~ra_iterator() = default;

    //! return element currently pointed at
    reference operator*() const
    {
        return host[pos];
    }

    //! comparison operators evaluate container values currently pointed at
    bool operator==(ra_iterator const & rhs) const
    {
        return (*this) == rhs;
    }

    bool operator!=(ra_iterator const & rhs) const
    {
        return (*this) != rhs;
    }

    bool operator<(ra_iterator const & rhs) const
    {
        return (*this) < rhs;
    }

    bool operator>(ra_iterator const & rhs) const
    {
        return (*this) > rhs;
    }

    bool operator<=(ra_iterator const & rhs) const
    {
        return (*this) <= rhs;
    }

    bool operator>=(ra_iterator const & rhs) const
    {
        return (*this) >= rhs;
    }

    //! prefix increment
    ra_iterator& operator++()
    {
        ++pos;
        return (*this);
    }

    //! postfix increment
    ra_iterator operator++(int)
    {
        ra_iterator ret{*this};
        ++(*this);
        return ret;
    }

    //! prefix decrement
    ra_iterator& operator--()
    {
        --pos;
        return *this;
    }

    //! postfix decrement
    ra_iterator operator--(int)
    {
        ra_iterator ret{*this};
        --(*this);
        return ret;
    }

    /*
    ra_iterator& operator+=(size_type); //optional
    ra_iterator operator+(size_type) const; //optional
    friend ra_iterator operator+(size_type, const ra_iterator&); //optional
    ra_iterator& operator-=(size_type); //optional
    ra_iterator operator-(size_type) const; //optional
    difference_type operator-(ra_iterator) const; //optional
    */

    pointer operator->() const
    {
        return &(*this);
    }

    // TODO: behaviour: with or without offset? range check here or delegation to container?
    reference operator[](size_type const n) const
    {
        return host[pos + n];
    }

};

template <typename container_type>
    requires container_concept<container_type> //random_access_range_concept<container_type>
struct ra_const_iterator
{
private:
    container_type const & host;
    typename container_type::size_type pos;
public:
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::value_type value_type;
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

    //! default constructor
    ra_const_iterator() : host(*(container_type*)0) {};
    //! initialize with pointer of container structure
    ra_const_iterator(container_type & _host) : host{_host}, pos{0} {};
    //! copy constructors
    ra_const_iterator(ra_const_iterator const &) = default;
    ra_const_iterator(ra_const_iterator &&) = default;

    ra_const_iterator& operator=(ra_const_iterator const &) = default;
    ~ra_const_iterator() = default;

    //! dereference operator
    // TODO: spec, range check?
    reference operator*() const
    {
        return host[pos];
    }

    bool operator==(ra_const_iterator const & rhs) const
    {
        return (*this) == rhs;
    }

    bool operator!=(ra_const_iterator const & rhs) const
    {
        return (*this) != rhs;
    }

    bool operator<(ra_const_iterator const & rhs) const
    {
        return (*this) < rhs;
    }

    bool operator>(ra_const_iterator const & rhs) const
    {
        return (*this) > rhs;
    }

    bool operator<=(ra_const_iterator const & rhs) const
    {
        return (*this) <= rhs;
    }

    bool operator>=(ra_const_iterator const & rhs) const
    {
        return (*this) >= rhs;
    }

    //! prefix increment
    ra_const_iterator& operator++()
    {
        ++pos;
        return (*this);
    }

    //! postfix increment
    ra_const_iterator operator++(int)
    {
        ra_iterator ret{*this};
        ++(*this);
        return ret;
    }

    //! prefix decrement
    ra_const_iterator& operator--()
    {
        --pos;
        return (*this);
    }

    //! postfix decrement
    ra_const_iterator operator--(int)
    {
        ra_iterator ret{*this};
        --(*this);
        return ret;
    }

    /*
    ra_const_iterator& operator+=(size_type); //optional
    ra_const_iterator operator+(size_type) const; //optional
    friend ra_const_iterator operator+(size_type, const ra_const_iterator&); //optional
    ra_const_iterator& operator-=(size_type); //optional
    ra_const_iterator operator-(size_type) const; //optional
    difference_type operator-(ra_const_iterator) const; //optional
    */

    pointer operator->() const
    {
        return &(*this);
    }

    // TODO: behaviour: with or without offset? range check here or delegation to container?
    reference operator[](size_type const n) const
    {
        return host[pos + n];
    }

};

} // namespace seqan3::detail
