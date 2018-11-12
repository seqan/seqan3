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
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief A constexpr string implementation to manipulate string literals at compile time.
 */

#pragma once

#include <array>
#include <type_traits>

namespace seqan3
{

/*!\brief Implements a constexpr string that can be used for compile time computations.
 * \ingroup container
 * \tparam N The size of the constexpr string.
 * \implements seqan3::ContainerRange
 * \implements std::ranges::RandomAccessRange
 *
 * This class provides a string type that can be constructed, evaluated and operated on at compile time.
 * The stored string can be accessed as either a std::string or a c-style string through the respective member functions.
 * All operations are constexpr except the conversion functions.
 *
 * Internally the string stores an array of size `N+1` containing the `0`-byte to allow conversion to
 * c-style strings. This also means that when creating an instance of this string from a string literal containing
 * a `0`-byte the size of this instance is one less.
 */
template <size_t N>
class constexpr_string
{
protected:
    //!\brief Alias for the underlying data type.
    using data_type = std::array<char, N + 1>;

    //!\brief Gives access to the merge constructor to constexpr_string with different size.
    template <size_t N2>
    friend class constexpr_string;

    //!\brief The internal string stored as array including \0-byte as last character.
    data_type lit;

    /*!\brief Constructs new constexpr_string by merging two other constexpr_strings.
     * \param lhs left-hand-side of constexpr_string to merge.
     * \param rhs right-hand-side of constexpr_string to merge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <size_t N1>
    constexpr constexpr_string(constexpr_string<N1> const & lhs,
                               constexpr_string<N - N1> const & rhs) noexcept
        : constexpr_string{}
    {
        for (unsigned i = 0; i < N1; ++i)
            lit[i] = lhs[i];
        for (unsigned i = N1; i < N + 1; ++i)
            lit[i] = rhs[i - N1];
    }

public:

    /*!\name Member types
     * \{
     */
    using value_type      = char;
    using reference       = char &;
    using const_reference = const char &;
    using iterator        = std::add_pointer_t<value_type>;
    using const_iterator  = std::add_pointer_t<value_type const>;
    using difference_type = typename data_type::difference_type;
    using size_type       = typename data_type::size_type;
    //!\}

    //!\cond
    // this signals to range-v3 that something is a container :|
    using allocator_type    = void;
    //!\endcond

    /*!\name Constructors, destructor and assignment
     * \{
     */
     //!\brief Default default constructor.
    constexpr constexpr_string() = default;
    //!\brief Default copy constructor.
    constexpr constexpr_string(constexpr_string const &) = default;
    //!\brief Default move constructor.
    constexpr constexpr_string(constexpr_string &&) = default;
    //!\brief Default copy-assignment operator.
    constexpr constexpr_string & operator=(constexpr_string const &) = default;
    //!\brief Default move-assignment operator.
    constexpr constexpr_string & operator=(constexpr_string &&) = default;

    /*!\brief Construction from literal.
     * \param _lit The literal to construct the string for.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr constexpr_string(const char (&_lit)[N + 1]) noexcept
        : constexpr_string{}
    {
        // static_assert(lit[N] == '\0'); TODO(rrahn): Fix me
        for (unsigned i = 0; i < N + 1; ++i)
            lit[i] = _lit[i];
    }

    /*!\brief Construction from char array.
     * \param src The char array to construct the constexpr_string for.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr constexpr_string(std::array<char, N> const & src) noexcept
        : constexpr_string{}
    {
        // static_assert(src[N-1] != '\0'); TODO(rrahn): Fix me
        for (unsigned i = 0; i < N; ++i)
            lit[i] = src[i];
        lit[N] = '\0';
    }

    /*!\brief Construction from char.
     * \param c The character to construct the constexpr_string for.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr constexpr_string(const char c) noexcept
        : constexpr_string{}
    {
        lit[0] = c;
        lit[1] = '\0';
    }

    //!\brief Default destructor.
    ~constexpr_string() = default;
    //!\}

    /*!\name Element access
     * \{
     */
     /*!\brief Access an element in the string.
      *
      * ### Exceptions
      *
      * No-throw guarantee.
      *
      * ### Complexity
      *
      * Constant.
      */
    constexpr reference operator[](size_t pos) noexcept
    {
        return lit[pos];
    }

    /*!\brief Access an element in the string.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Constant.
     */
    constexpr const_reference operator[](size_t pos) const noexcept
    {
        return lit[pos];
    }

    /*!\brief Returns the content represented as std::string.
     *
     * \returns `std::string` The stored string.
     *
     * ### Exceptions
     *
     * Strong exception guarantee. No data is modified.
     *
     * ### Complexity
     *
     * Linear in the size of the string.
     */
    std::string string() const
    {
        return std::string{ lit.cbegin(), lit.cend() - 1};
    }

    /*!\brief Returns the content represented as 0-terminated c-style string.
     *
     * \returns `const char *` The stored string.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Constant.
     */
    constexpr const char * c_str() const noexcept
    {
        return lit.data();
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    //!\brief Returns the size.
    constexpr size_type size() const noexcept
    {
        return N;
    }

    //!\brief Returns the maximal capacity (same as size()).
    constexpr size_type max_size() const noexcept
    {
        return size();
    }

    //!\brief Determines whether the string is empty.
    constexpr bool empty() const noexcept
    {
        return cbegin() == cend();
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the begin to the string.
    constexpr iterator begin() noexcept
    {
        return &lit[0];
    }

    //!\copydoc constexpr_string::begin().
    constexpr const_iterator begin() const noexcept
    {
        return &lit[0];
    }

    //!\copydoc constexpr_string::begin().
    constexpr const_iterator cbegin() const noexcept
    {
        return &lit[0];
    }

    //!\brief Returns iterator pass the end of the string.
    constexpr iterator end() noexcept
    {
        return &lit[N];
    }

    //!\copydoc constexpr_string::end().
    constexpr const_iterator end() const noexcept
    {
        return &lit[N];
    }

    //!\copydoc constexpr_string::end().
    constexpr const_iterator cend() const noexcept
    {
        return &lit[N];
    }
    //!\}

    /*!\name Operations
     * \{
     */

    /*!\brief Concatenates two constexpr_strings.
     * \param rhs The right-hand-side to concat with.
     * \returns `constexpr_string<N + N2>` The new constexpr_string with size `N` + `N2`.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Linear in the size of the strings.
     */
    template <size_t N2>
    constexpr constexpr_string<N + N2> operator+(constexpr_string<N2> const & rhs) const noexcept
    {
        return constexpr_string<N + N2>{*this, rhs};
    }

    /*!\brief Swaps the contents.
     *
     * Both strings must have the same size in order to swap them.
     *
     * ### Complexity
     *
     * Constant.
     */
    constexpr void swap(constexpr_string & other) noexcept
    {
        std::swap(*this, other);
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    /*!\brief Compares two strings lexicographically.
     * \returns `bool` `true` if the corresponding comparison holds, `false` otherwise.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Linear in the size of the strings.
     */
    template <size_t N2>
    constexpr bool operator==(constexpr_string<N2> const & rhs) const noexcept
    {
        if constexpr (N != N2)
            return false;

        const_iterator it_rhs = rhs.cbegin();
        for (const_iterator it = cbegin(); it != cend(); ++it, ++it_rhs)
        {
            if (*it != *it_rhs)
                return false;
        }

        return true;
    }

    //!\copydoc constexpr_string::operator==()
    template <size_t N2>
    constexpr bool operator!=(constexpr_string<N2> const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    //!\copydoc constexpr_string::operator==()
    template <size_t N2>
    constexpr bool operator<(constexpr_string<N2> const & rhs) const noexcept
    {
        for (unsigned i = 0; i < ((N < N2) ? N : N2); ++i)
        {
            if (lit[i] < rhs.lit[i])
                return true;
            else if (lit[i] != rhs.lit[i])
                return false;
        }
        return N < N2;
    }

    //!\copydoc constexpr_string::operator==()
    template <size_t N2>
    constexpr bool operator<=(constexpr_string<N2> const & rhs) const noexcept
    {
        for (unsigned i = 0; i < ((N < N2) ? N : N2); ++i)
        {
            if (lit[i] > rhs.lit[i])
                return false;
        }
        return N <= N2;
    }

    //!\copydoc constexpr_string::operator==()
    template <size_t N2>
    constexpr bool operator>(constexpr_string<N2> const & rhs) const noexcept
    {
        return !(*this <= rhs);
    }

    //!\copydoc constexpr_string::operator==()
    template <size_t N2>
    constexpr bool operator>=(constexpr_string<N2> const & rhs) const noexcept
    {
        return !(*this < rhs);
    }
    //!\}
};

/*!\name Operations
 * \{
 */
/*!\brief Exchanges the given values.
 * \relates constexpr_string
 * \tparam N The size of the constexpr_string.
 * \param lhs The left operand.
 * \param rhs The right operand.
 *
 * ### Exceptions
 *
 * No-throw guarantee.
 *
 * ### Complexity
 *
 * Constant.
 */
template <size_t N>
constexpr inline void
swap(constexpr_string<N> & lhs, constexpr_string<N> & rhs)
{
    lhs.swap(rhs);
}
//!\}

/*!\name Deduction guides
 * \{
 */
//!\brief Deduces constexpr_string from string literals.
//!\relates constexpr_string
template <size_t N>
constexpr_string(const char (&)[N]) -> constexpr_string<N - 1>;

//!\brief Deduces constexpr_string from std::array of type char.
//!\relates constexpr_string
template <size_t N>
constexpr_string(std::array<char, N> const &) -> constexpr_string<N>;

//!\brief Deduces constexpr_string from char.
//!\relates constexpr_string
constexpr_string(const char) -> constexpr_string<1>;
//!\}
}  // namespace seqan3
