// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief A constexpr string implementation to manipulate string literals at compile time.
 */

#pragma once

#include <seqan3/range/container/small_vector.hpp>

namespace seqan3
{

/*!\brief Implements a small string that can be used for compile time computations.
 * \ingroup container
 * \implements seqan3::reservible_container
 * \tparam capacity_ The capacity of the small string.
 *
 * This class provides a string type for small strings and compile-time contexts. It has fixed capacity, but variable
 * size within the capacity. It is always allocated on the stack and most of it's members are `constexpr`-qualified.
 * The underlying data can be exposed as a null-terminated c-style string (without copying) and conversion operators to
 *  std::string are provided (this involves copying).
 *
 * ### Implementation notes
 *
 * Internally the string stores a null-terminated array of size `capacity_ + 1` and the size of the string as a member.
 * The smallest possible type is used for storage of the size. For example, `small_string<30>` uses 32bytes of memory
 * (one byte extra for the null-terminator and one byte to save the size).
 *
 * Usage:
 *
 * \include test/snippet/range/container/small_string.cpp
 *
 */
template <size_t capacity_>
class small_string : public small_vector<char, capacity_ + 1>
{
private:
    //!\brief The parent class.
    using base_t = small_vector<char, capacity_ + 1>;

    // make data inherited members visible
    using base_t::data_;
    using base_t::sz;
public:
    /*!\name Associated types
     * \{
     */
    using typename base_t::value_type;
    using typename base_t::reference;
    using typename base_t::const_reference;
    using typename base_t::iterator;
    using typename base_t::const_iterator;
    using typename base_t::difference_type;
    using typename base_t::size_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Inherit the parent class constructors.
    using base_t::base_t;
    //!\brief Inherit the parent class member function assign.
    using base_t::assign;

    /*!\brief Construction from literal.
     * \param _lit The literal to construct the string for.
     *
     * The `char` literal is expected to be null-terminated (asserted in debug-mode). If it is not, the last character
     * will be lost when copying to the instance of `small_string`.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <size_t N>
    constexpr small_string(char const (&_lit)[N]) noexcept : small_string{}
    {
        static_assert(N <= capacity_ + 1, "Length of string literal exceeds capacity of small_string.");
        assign(_lit);
    }

    /*!\brief Construction from char.
     * \param c The character to construct the small_string for.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    explicit constexpr small_string(char const c) noexcept : small_string{}
    {
        assign(1, c);
    }

    /*!\brief Assign from literal.
     * \param _lit The literal to assign the string from.
     *
     * The `char` literal is expected to be null-terminated (asserted in debug-mode). If it is not, the last character
     * will be lost when copying to the instance of `small_string`.
     *
     * ### Complexity
     *
     * Linear in the size of _lit.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <size_t N>
    constexpr small_string & operator=(char const (&_lit)[N]) noexcept
    {
        static_assert(N <= capacity_ + 1, "Length of string literal exceeds capacity of small_string.");
        assign(_lit);
        return *this;
    }

    /*!\brief Assign from literal.
     * \param _lit The literal to assign the string from.
     *
     * The `char` literal is expected to be null-terminated (asserted in debug-mode). If it is not, the last character
     * will be lost when copying to the instance of `small_string`.
     *
     * ### Complexity
     *
     * Linear in the size of _lit.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <size_t N>
    constexpr void assign(char const (&_lit)[N]) noexcept
    {
        static_assert(N <= capacity_ + 1, "Length of string literal exceeds capacity of small_string.");
        assert(_lit[N - 1] == '\0');
        assign(&_lit[0], &_lit[N-1]);
    }

    /*!\brief Assign from pair of iterators.
     * \tparam begin_it_type Must satisfy std::forward_iterator and the `value_type` must be constructible from
     *                       the reference type of begin_it_type.
     * \tparam   end_it_type Must satisfy std::sentinel_for.
     * \param[in]   begin_it Begin of range to construct/assign from.
     * \param[in]     end_it End of range to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the distance between `begin_it` and `end_it`.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     */
    template <std::forward_iterator begin_it_type, typename end_it_type>
    //!\cond
        requires std::sentinel_for<end_it_type, begin_it_type> &&
                 std::constructible_from<value_type, std::iter_reference_t<begin_it_type>>
    //!\endcond
    constexpr void assign(begin_it_type begin_it, end_it_type end_it) noexcept
    {
        base_t::assign(begin_it, end_it);
        data_[sz] = '\0';
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    //!\brief Returns the maximal size which equals the capacity.
    static constexpr size_type max_size() noexcept
    {
        return capacity_;
    }

    //!\brief Returns the maximal capacity.
    static constexpr size_type capacity() noexcept
    {
        return capacity_;
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    //!\copydoc seqan3::small_vector::clear
    constexpr void clear() noexcept
    {
        sz = 0;
        data_[0] = '\0';
    }

    //!\copydoc seqan3::small_vector::push_back(value_type const value)
    constexpr void push_back(char const value) noexcept
    {
        assert(sz < capacity_);
        data_[sz] = value;
        ++sz;
        data_[sz] = '\0';
    }

    //!\copydoc seqan3::small_vector::pop_back
    constexpr void pop_back() noexcept
    {
        assert(sz > 0);
        --sz;
        data_[sz] = '\0';
    }

    //!\copydoc seqan3::small_vector::resize(size_type const)
    constexpr void resize(size_type const count) noexcept
    {
        resize(count, '\0');
    }

    //!\copydoc seqan3::small_vector::resize(size_type const, value_type const value)
    constexpr void resize(size_type const count, char const value) noexcept
    {
        assert(count <= capacity_);

        for (size_t i = sz; i < count; ++i)  // sz < count; add `value` in [sz, count)
            data_[i] = value;

        sz = count;
        data_[sz] = '\0';
    }

    /*!\brief Removes specified elements from the container.
     * \param   index Remove the elements starting at `index`. Defaults to `0`.
     * \param   count The number of elements to remove. Defaults to `max_size()`.
     * \returns `*this`
     *
     * Invalidates iterators and references at or after the point of the erase, including the end() iterator.
     *
     * The iterator `pos` must be valid and dereferenceable. Thus the end() iterator (which is valid, but is not
     * dereferencable) cannot be used as a value for pos.
     *
     * ### Complexity
     *
     * Linear in size().
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr small_string & erase(size_type index = 0, size_type count = max_size()) noexcept
    {
        assert(index <= this->size());

        iterator it = this->begin() + index;
        base_t::erase(it, it + std::min<size_type>(count, this->size() - index));
        return *this;
    }
    //!\}

    /*!\name Concatenation
     * \{
     */
    /*!\brief Concatenates two small_strings by returning a new small_string.
     * \param lhs The left-hand-side to concat with.
     * \param rhs The right-hand-side to concat with.
     * \returns `small_string<capacity_ + capacity2>` The new small_string with size `capacity_` + `capacity2`.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Linear in the size of the strings.
     */
    template <size_t capacity2>
    constexpr friend small_string<capacity_ + capacity2> operator+(small_string const & lhs,
                                                                   small_string<capacity2> const & rhs) noexcept
    {
        small_string<capacity_ + capacity2> tmp{lhs};
        tmp.insert(tmp.end(), rhs.begin(), rhs.end());
        return tmp;
    }
    //!\}

    /*!\name Conversion
     * \{
     */

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
    std::string str() const
    {
        return std::string{this->cbegin(), this->cend()};
    }

    /*!\brief Returns the content represented as 0-terminated c-style string.
     *
     * \returns `char const *` The stored string.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Constant.
     */
    constexpr char const * c_str() const noexcept
    {
        return data_.data();
    }

    /*!\brief Implicit conversion to std::string which delegates to seqan3::small_string::str().
     *
     * ### Exceptions
     *
     * Strong exception guarantee. No data is modified.
     *
     * ### Complexity
     *
     * Linear in the size of the string.
     */
    operator std::string() const
    {
        return str();
    }
    //!\}

    /*!\name Input/output
     * \{
     */

    /*!\brief Formatted output for the seqan3::small_string.
     * \param[in,out] os  The std::basic_ostream to write to.
     * \param[in]     str The seqan3::small_string to read from.
     * \returns `os`.
     *
     * \details
     *
     * Internally calls `os << str.str()`.
     */
    friend std::ostream & operator<<(std::ostream & os, small_string const & str)
    {
        os << str.str();
        return os;
    }

    /*!\brief Formatted input for the seqan3::small_string.
     * \param[in,out] is  The std::basic_istream to read from.
     * \param[out]    str The seqan3::small_string to write to.
     * \returns `is`.
     *
     * \details
     *
     * Reads at most seqan3::small_string::max_size characters from the stream.
     * If a stream error occurred or no characters could be extracted the std::ios_base::failbit is set.
     * This may throw an exception.
     */
    friend std::istream & operator>>(std::istream & is, small_string & str)
    {
        // Check if stream is ok and skip leading whitespaces.
        std::istream::sentry s(is);
        if (s)
        {
            str.erase(); // clear the string
            std::streamsize num_char = (is.width() > 0)
                ? std::min<std::streamsize>(is.width(), str.max_size())
                : str.max_size();
            assert(num_char > 0);
            for (std::streamsize n = num_char; n > 0 && !std::isspace(static_cast<char>(is.peek()), is.getloc()); --n)
            {
                char c = is.get();
                if (is.eof())
                    break;
                str.push_back(c);
            }

            if (str.size() == 0) // nothing extracted so we set the fail bit.
                is.setstate(std::ios_base::failbit);

            is.width(0); // cancel the effects of std::setw, if any.
        }

        return is;
    }
    //!\}
};

/*!\name Deduction guides
 * \{
 */
//!\brief Deduces small_string from string literals.
//!\relates small_string
template <size_t N>
small_string(char const (&)[N]) -> small_string<N - 1>;

//!\brief Deduces small_string from std::array of type char.
//!\relates small_string
template <size_t N>
small_string(std::array<char, N> const &) -> small_string<N>;

//!\brief Deduces small_string from char.
//!\relates small_string
small_string(char const) -> small_string<1>;
//!\}

} // namespace seqan3
