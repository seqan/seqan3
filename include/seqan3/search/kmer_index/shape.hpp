// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \author Vinzenz May <vinzenz.may AT fu-berlin.de>
 * \brief Provides seqan3::shape.
 */

#pragma once

#include <seqan3/utility/container/dynamic_bitset.hpp>

namespace seqan3
{

/*!\brief A strong type of underlying type `uint8_t` that represents the ungapped shape size.
 * \ingroup search_kmer_index
 */
struct ungapped
{
    //!\brief The ungapped shape size.
    uint8_t value;
};

/*!\brief A strong type of underlying type `uint64_t` that represents the shape in binary representation.
 * \ingroup search_kmer_index
 */
struct bin_literal
{
    //!\brief The shape in binary representation.
    uint64_t value;
};

/*!\brief A class that defines which positions of a pattern to hash.
 * \ingroup search_kmer_index
 *
 * \details
 *
 * \experimentalapi
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
    constexpr shape() noexcept = default;                          //!< Defaulted.
    constexpr shape(shape const &) noexcept = default;             //!< Defaulted.
    constexpr shape(shape &&) noexcept = default;                  //!< Defaulted.
    constexpr shape & operator=(shape const &) noexcept = default; //!< Defaulted.
    constexpr shape & operator=(shape &&) noexcept = default;      //!< Defaulted.
    ~shape() noexcept = default;                                   //!< Defaulted.

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

inline namespace literals
{

/*!\name Other literals
 * \{
 */
/*!\brief The seqan3::shape literal.
 * \param[in] value The unsigned integer to assign.
 * \relatesalso seqan3::shape
 * \returns seqan3::shape
 */
constexpr shape operator""_shape(unsigned long long const value)
{
    return shape{bin_literal{value}};
}
//!\}

} // namespace literals

} // namespace seqan3
