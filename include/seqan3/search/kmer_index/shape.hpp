// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \author Vinzenz May <vinzenz.may AT fu-berlin.de>
 * \brief Provides seqan3::shape.
 */

#pragma once

#include <array>

#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief A class that defines what positions of a pattern to hash.
*
* A seqan3::shape...
* \todo snippet
*
* \attention 0 < size <= 32
*/

class shape : public std::array<bool, 32>
{
private:
    //!\brief The size of the shape.
    size_t k;

public:

    using std::array<bool, 32>::array;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    shape()                                    noexcept = default; //!< Defaulted.
    constexpr shape(shape const &)             noexcept = default; //!< Defaulted.
    constexpr shape(shape &&)                  noexcept = default; //!< Defaulted.
    constexpr shape & operator=(shape const &) noexcept = default; //!< Defaulted.
    constexpr shape & operator=(shape &&)      noexcept = default; //!< Defaulted.
    ~shape()                                   noexcept = default; //!< Defaulted.

    /*!\brief Construct from a given size.
    * ### Complexity
    *
    * Constant.
    *
    * ### Exceptions
    *
    * No-throw guarantee.
    */
    constexpr shape(size_t k_) noexcept : k(k_)
    {
        assert(k > 0 && k < 33);
        fill(false);
        std::fill_n(begin(), k_, true);
    }

    /*!\brief Construct from a parameter pack.
    * ### Complexity
    *
    * Constant
    *
    * ### Exceptions
    *
    * No-throw guarantee.
    */
    template <typename... t>
    //!\cond
        requires sizeof...(t) >= 2 && sizeof...(t) <= 32 && (std::Integral<t> && ...)
    //!\endcond
    constexpr shape(t... args) noexcept : k(sizeof...(args))
    {
        size_t i = 0;
        for (bool a : { args... })
            (*this)[i++] = a;
    }
    //!\}

    //!\brief The size of the shape.
    constexpr size_t size() const noexcept
    {
        return k;
    }
};

}// namespace seqan3
