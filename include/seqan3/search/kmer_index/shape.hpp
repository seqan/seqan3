// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \author Vinzenz May <vinzenz.may AT fu-berlin.de>
 * \brief Provides seqan3::shape.
 */

#pragma once

#include <seqan3/range/container/dynamic_bitset.hpp>

namespace seqan3
{

/*!\brief A strong type of underlying type `uint8_t` that represents the ungapped shape size.
 * \ingroup submodule_kmer_index
 */
struct ungapped
{
    //!\brief The ungapped shape size.
    uint8_t value;
};

/*!\brief A strong type of underlying type `uint64_t` that represents the shape in binary representation.
 * \ingroup submodule_kmer_index
 */
struct bin_literal
{
    //!\brief The shape in binary representation.
    uint64_t value;
};

/*!\brief A class that defines what positions of a pattern to hash.
 * \ingroup submodule_kmer_index
 *
 * \details
 *
 * When hashing a sequence, there may be positions that do not count towards the final hash value.
 * A shape offers an easy way to define such patterns. Given a k-mer length `k` (0 < `k` <= 58), a shape
 * represents a binary sequence where a `0` encodes a "don't care position", i.e. a position that is not taken into
 * account when computing the hash value. A `1` therefore translates to a position that is used to compute the hash
 * value.
 *
 * ### Example
 *
 * \include test/snippet/search/kmer_index/shape.cpp
 *
 * \attention 0 < size <= 58
 */

class shape : public dynamic_bitset<58>
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr shape()                          noexcept = default; //!< Defaulted.
    constexpr shape(shape const &)             noexcept = default; //!< Defaulted.
    constexpr shape(shape &&)                  noexcept = default; //!< Defaulted.
    constexpr shape & operator=(shape const &) noexcept = default; //!< Defaulted.
    constexpr shape & operator=(shape &&)      noexcept = default; //!< Defaulted.
    ~shape()                                   noexcept = default; //!< Defaulted.

    /*!\brief Construct an ungapped shape from a given size.
     *
     * \details
     *
     * ### Complexity
     *
     * Linear in `k`.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * \attention The size must be in the interval [1, 58]. In Debug mode, an assertion will check this constraint.
     */
    constexpr shape(ungapped k) noexcept : dynamic_bitset<58>((1ULL << k.value) - 1)
    {
        assert(k.value > 0);
    }

    /*!\brief Construct from a given seqan3::bin_literal.
     *
     * \details
     *
     * ### Complexity
     *
     * Linear in the size of the `bin_literal`.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * \attention The size of the `bin_literal` must be in the interval [1, 58].
     */
    constexpr shape(bin_literal const literal) noexcept : dynamic_bitset<58>(literal.value)
    {
        assert(front() == 1); // First position must be 1, e.g. no 0111 shape
        assert(back() == 1);  // Last position must be 1, e.g. no 1110 shape
    }
    //!\}
};

/*!\brief The seqan3::shape literal.
 * \param[in] value The unsigned integer to assign.
 * \relates seqan3::shape
 * \returns seqan3::shape
 */
constexpr shape operator""_shape(unsigned long long const value)
{
    return shape{bin_literal{value}};
}

}// namespace seqan3
