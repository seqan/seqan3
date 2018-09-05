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
 * \brief Provides seqan3::bitcompressed_vector.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <range/v3/iterator_range.hpp>

#include <sdsl/int_vector.hpp>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief A space-optimised version of std::vector that compresses multiple letters into a single byte.
 * \tparam alphabet_type The value type of the container, must satisfy seqan3::alphabet_concept and not be `&`.
 * \implements seqan3::reservable_container_concept
 * \implements seqan3::cerealisable_concept
 * \ingroup container
 *
 * This class template behaves just like std::vector<alphabet_type> but has an internal representation where
 * multiple values are packed into a single byte/word to save space, e.g. bitcompressed_vector<seqan3::dna4> uses a
 * quarter of of the memory that std::vector<seqan3::dna4> uses, because a single seqan3::dna4 letter can be represented
 * in two bits (instead of 8 which is the lower bound for a single object in C++).
 *
 * The disadvantages are slightly slower operations and unsafety towards parallel writes to adjacent positions
 * in the bitcompressed_vector.
 *
 * ### Example
 *
 * \snippet test/snippet/range/container/bitcompressed_vector.cpp usage
 *
 * ### Thread safety
 *
 * This container provides no thread-safety beyond the promise given also by the STL that all
 * calls to `const` member function are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * An important difference to std::vector is that calling `vec[i] = value;` and `vec[j] = value2;` from two different
 * threads at the same time **is not safe** and will lead to corruption if both values are stored in the same
 * 64bit-block, i.e. if the distance between `i` and `j` is smaller than 64 / alphabet_size.
 */
template <alphabet_concept alphabet_type>
//!\cond
    requires std::is_same_v<alphabet_type, std::remove_reference_t<alphabet_type>>
//!\endcond
class bitcompressed_vector
{
private:
    //!\brief The number of bits needed to represent a single letter of the alphabet_type.
    static constexpr size_t bits_per_letter = std::ceil(std::log2(alphabet_size_v<alphabet_type>));

    static_assert(bits_per_letter <= 64, "Alphabet must be representable in at most 64bit.");

    //!\brief Type of the underlying SDSL vector.
    using data_type = sdsl::int_vector<bits_per_letter>;

    //!\brief The data storage.
    data_type data;

    //!\brief Proxy data type returned by seqan3::bitcompressed_vector as reference to element.
    class reference_proxy_type
    {
    private:
        //!\brief The proxy of the underlying data type.
        sdsl::int_vector_reference<data_type> internal_proxy;
    public:
        /*!\name Constructors, destructor and assignment
         * \brief All are explicitly defaulted.
         * \{
         */
        reference_proxy_type() = default;
        constexpr reference_proxy_type(reference_proxy_type const &) = default;
        constexpr reference_proxy_type(reference_proxy_type &&) = default;
        constexpr reference_proxy_type & operator=(reference_proxy_type const &) = default;
        constexpr reference_proxy_type & operator=(reference_proxy_type &&) = default;
        ~reference_proxy_type() = default;

        //!\brief Initialise from internal proxy type.
        reference_proxy_type(sdsl::int_vector_reference<data_type> const & internal) :
            internal_proxy{internal}
        {}
        //!\}

        /*!\name Conversion from/to alphabet_type
         * \{
         */
        //!\brief Implicitly convertible to alphabet_type.
        constexpr operator alphabet_type() const
            noexcept(noexcept(assign_rank(alphabet_type{}, static_cast<uint64_t>(internal_proxy))))
        {
            return assign_rank(alphabet_type{}, static_cast<uint64_t>(internal_proxy));
        }

        //!\brief Assignable from alphabet_type.
        constexpr reference_proxy_type & operator=(alphabet_type const a)
            noexcept(noexcept(internal_proxy = to_rank(a)))
        {
            internal_proxy = to_rank(a);
            return *this;
        }
        //!\}

        //!\name Comparison operators (self)
        //!\{
        constexpr bool operator==(reference_proxy_type const & rhs) const noexcept
        {
            return internal_proxy == rhs.internal_proxy;
        }

        constexpr bool operator!=(reference_proxy_type const & rhs) const noexcept
        {
            return !(internal_proxy == rhs.internal_proxy);
        }

        constexpr bool operator<(reference_proxy_type const & rhs) const noexcept
        {
            return internal_proxy < rhs.internal_proxy;
        }

        constexpr bool operator>(reference_proxy_type const & rhs) const noexcept
        {
            return !(internal_proxy <= rhs.internal_proxy);
        }

        constexpr bool operator<=(reference_proxy_type const & rhs) const noexcept
        {
            return internal_proxy == rhs.internal_proxy || internal_proxy < rhs.internal_proxy;
        }

        constexpr bool operator>=(reference_proxy_type const & rhs) const noexcept
        {
            return !(internal_proxy < rhs.internal_proxy);
        }
        //!\}

        //!\name Comparison operators (alphabet_type)
        //!\{
        constexpr bool operator==(alphabet_type const & rhs) const noexcept
        {
            return static_cast<alphabet_type>(*this) == rhs;
        }

        constexpr bool operator!=(alphabet_type const & rhs) const noexcept
        {
            return !(static_cast<alphabet_type>(*this) == rhs);
        }

        constexpr bool operator<(alphabet_type const & rhs) const noexcept
        {
            return static_cast<alphabet_type>(*this) < rhs;
        }

        constexpr bool operator>(alphabet_type const & rhs) const noexcept
        {
            return !(static_cast<alphabet_type>(*this) <= rhs);
        }

        constexpr bool operator<=(alphabet_type const & rhs) const noexcept
        {
            return static_cast<alphabet_type>(*this) == rhs || static_cast<alphabet_type>(*this) < rhs;
        }

        constexpr bool operator>=(alphabet_type const & rhs) const noexcept
        {
            return !(static_cast<alphabet_type>(*this) < rhs);
        }
        //!\}

        static_assert(!std::WeaklyIncrementable<alphabet_type>,
                      "TODO: bitcompressed_vector::reference_proxy_type needs to be adapted to also be incrementable.");
    };

    //!\cond
    //NOTE(h-2): it is entirely unclear to me why we need this
    template <typename t>
        requires std::is_same_v<value_type_t<remove_cvref_t<t>>, alphabet_type>
    static constexpr bool has_same_value_type_v = true;
    //!\endcond

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Equals the alphabet_type.
    using value_type        = alphabet_type;
    //!\brief A proxy type that enables assignment, if the underlying data structure also provides a proxy.
    using reference         = std::conditional_t<std::is_lvalue_reference_v<reference_t<data_type>>,
                                                 reference_t<data_type>,
                                                 reference_proxy_type>;
    //!\brief Equals the alphabet_type / value_type.
    using const_reference   = alphabet_type;
    //!\brief The iterator type of this container (a random access iterator).
    using iterator          = detail::random_access_iterator<bitcompressed_vector>;
    //!\brief The const_iterator type of this container (a random access iterator).
    using const_iterator    = detail::random_access_iterator<bitcompressed_vector const>;
    //!\brief A signed integer type (usually std::ptrdiff_t)
    using difference_type   = ranges::difference_type_t<data_type>;
    //!\brief An unsigned integer type (usually std::size_t)
    using size_type         = ranges::size_type_t<data_type>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bitcompressed_vector() = default;
    constexpr bitcompressed_vector(bitcompressed_vector const &) = default;
    constexpr bitcompressed_vector(bitcompressed_vector &&) = default;
    constexpr bitcompressed_vector & operator=(bitcompressed_vector const &) = default;
    constexpr bitcompressed_vector & operator=(bitcompressed_vector &&) = default;
    ~bitcompressed_vector() = default;

    /*!\brief Construct from a different range.
     * \tparam other_range_t The type of range to construct from; must satisfy std::ranges::InputRange and
     *                       std::CommonReference<value_type_t<other_range_t>, value_type>.
     * \param[in]      range The sequences to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the size of `range`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::InputRange other_range_t>
    //!\cond
        requires has_same_value_type_v<other_range_t>
    //!\endcond
    explicit bitcompressed_vector(other_range_t && range) :
        bitcompressed_vector{ranges::begin(range), ranges::end(range)}
    {}

    /*!\brief Construct with `count` times `value`.
     * \param[in] count Number of elements.
     * \param[in] value The initial value to be assigned.
     *
     * ### Complexity
     *
     * In \f$O(count)\f$.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    bitcompressed_vector(size_type const count, value_type const value) :
        data(count, to_rank(value))
    {}

    /*!\brief Construct from pair of iterators.
     * \tparam begin_iterator_type Must model std::ForwardIterator and
     *                             std::CommonReference<value_type_t<begin_iterator_type>, value_type>.
     * \tparam   end_iterator_type Must model std::Sentinel.
     * \param[in]         begin_it Begin of range to construct/assign from.
     * \param[in]           end_it End of range to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the distance between `begin_it` and `end_it`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ForwardIterator begin_iterator_type, std::Sentinel<begin_iterator_type> end_iterator_type>
    bitcompressed_vector(begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires std::CommonReference<value_type_t<begin_iterator_type>, value_type>
    //!\endcond
    {
        insert(cend(), begin_it, end_it);
    }

    /*!\brief Construct from `std::initializer_list`.
     * \param[in] ilist A `std::initializer_list` of value_type.
     *
     * ### Complexity
     *
     * Linear in the size of `ilist`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    bitcompressed_vector(std::initializer_list<value_type> ilist) :
        bitcompressed_vector(std::begin(ilist), std::end(ilist))
    {}

    /*!\brief Assign from `std::initializer_list`.
     * \param[in] ilist A `std::initializer_list` of value_type.
     *
     * ### Complexity
     *
     * Linear in the size of `ilist`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    bitcompressed_vector & operator=(std::initializer_list<value_type> ilist)
    {
        assign(std::begin(ilist), std::end(ilist));
        return *this;
    }

    /*!\brief Assign from a different range.
     * \tparam other_range_t The type of range to be inserted; must satisfy std::ranges::InputRange and
     *                       std::CommonReference<value_type_t<other_range_t>, value_type>.
     * \param[in]      range The sequences to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the size of `range`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::InputRange other_range_t>
    void assign(other_range_t && range)
    //!\cond
        requires std::CommonReference<value_type_t<other_range_t>, value_type>
    //!\endcond
    {
        bitcompressed_vector rhs{std::forward<other_range_t>(range)};
        swap(rhs);
    }

    /*!\brief Assign with `count` times `value`.
     * \param[in] count Number of elements.
     * \param[in] value The initial value to be assigned.
     *
     * ### Complexity
     *
     * In \f$O(count)\f$.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void assign(size_type const count, value_type const value)
    {
        bitcompressed_vector rhs{count, value};
        swap(rhs);
    }

    /*!\brief Assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy std::ForwardIterator and
     *                             std::CommonReference<value_type_t<begin_iterator_type>, value_type>.
     * \tparam   end_iterator_type Must satisfy std::Sentinel.
     * \param[in]         begin_it Begin of range to construct/assign from.
     * \param[in]           end_it End of range to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the distance between `begin_it` and `end_it`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ForwardIterator begin_iterator_type, std::Sentinel<begin_iterator_type> end_iterator_type>
    void assign(begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires std::CommonReference<value_type_t<begin_iterator_type>, value_type>
    //!\endcond
    {
        bitcompressed_vector rhs{begin_it, end_it};
        swap(rhs);
    }

    /*!\brief Assign from `std::initializer_list`.
     * \param[in] ilist A `std::initializer_list` of value_type.
     *
     * ### Complexity
     *
     * Linear in the size of `ilist`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void assign(std::initializer_list<value_type> ilist)
    {
        assign(std::begin(ilist), std::end(ilist));
    }

    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return iterator{*this};
    }

    //!\copydoc begin()
    const_iterator begin() const noexcept
    {
        return const_iterator{*this};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this};
    }

    /*!\brief Returns an iterator to the element following the last element of the container.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator end() noexcept
    {
        return iterator{*this, size()};
    }

    //!\copydoc end()
    const_iterator end() const noexcept
    {
        return const_iterator{*this, size()};
    }

    //!\copydoc end()
    const_iterator cend() const noexcept
    {
        return const_iterator{*this, size()};
    }
    //!\}

    /*!\name Element access
     * \{
     */
    /*!\brief Return the i-th element.
     * \param[in]              i Index of the element to retrieve.
     * \throws std::out_of_range If you access an element behind the last.
     * \returns                  Either a writable proxy to the element or a copy (if called in const context).
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws std::out_of_range if `i >= size()`.
     */
    reference at(size_type const i)
    {
        if (i >= size()) // [[unlikely]]
        {
            throw std::out_of_range{"Trying to access element behind the last in bitcompressed_vector."};
        }
        return (*this)[i];
    }

    //!\copydoc at()
    const_reference at(size_type const i) const
    {
        if (i >= size()) // [[unlikely]]
        {
            throw std::out_of_range{"Trying to access element behind the last in bitcompressed_vector."};
        }
        return (*this)[i];
    }

    /*!\brief Return the i-th element.
     * \param i The element to retrieve.
     * \returns Either a writable proxy to the element or a copy (if called in const context).
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    reference operator[](size_type const i) noexcept
    {
        assert(i < size());
        return data[i];
    }

    //!\copydoc operator[]()
    const_reference operator[](size_type const i) const noexcept
    {
        assert(i < size());
        return assign_rank(const_reference{}, data[i]);
    }

    /*!\brief Return the first element. Calling front on an empty container is undefined.
     * \returns Either a writable proxy to the element or a copy (if called in const context).
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    reference front() noexcept
    {
        assert(size() > 0);
        return (*this)[0];
    }

    //!\copydoc front()
    const_reference front() const noexcept
    {
        assert(size() > 0);
        return (*this)[0];
    }

    /*!\brief Return the last element.
     * \returns Either a writable proxy to the element or a copy (if called in const context).
     *
     * Calling back on an empty container is undefined. In debug mode an assertion checks the size of the container.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    reference back() noexcept
    {
        assert(size() > 0);
        return (*this)[size()-1];
    }

    //!\copydoc back()
    const_reference back() const noexcept
    {
        assert(size() > 0);
        return (*this)[size()-1];
    }

    /*!\name Capacity
     * \{
     */
    /*!\brief Checks whether the container is empty.
     * \returns `true` if the container is empty, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool empty() const noexcept
    {
        return size() == 0;
    }

    /*!\brief Returns the number of elements in the container, i.e. std::distance(begin(), end()).
     * \returns The number of elements in the container.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
    {
        return data.size();
    }

    /*!\brief Returns the maximum number of elements the container is able to hold due to system or library
     * implementation limitations, i.e. std::distance(begin(), end()) for the largest container.
     * \returns The number of elements in the container.
     *
     * This value typically reflects the theoretical limit on the size of the container. At runtime, the size
     * of the container may be limited to a value smaller than max_size() by the amount of RAM available.
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type max_size() const noexcept
    {
        return data.max_size();
    }

    /*!\brief Returns the number of elements that the container has currently allocated space for.
     * \returns The capacity of the currently allocated storage.
     *
     * \attention
     *
     * This does not operate on underlying concat container, see concat_capacity().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type capacity() const noexcept
    {
        return data.capacity();
    }

    /*!\brief Increase the capacity to a value that's greater or equal to new_cap.
     * \param            new_cap The new capacity.
     * \throws std::length_error If new_cap > max_size().
     * \throws    std::exception Any exception thrown by `Allocator::allocate()` (typically `std::bad_alloc`).
     *
     * Increase the capacity of the vector to a value that's greater or equal to new_cap.
     * If new_cap is greater than the current capacity(), new storage is allocated, otherwise the method does nothing.
     * If new_cap is greater than capacity(), all iterators, including the past-the-end iterator, and all references
     * to the elements are invalidated. Otherwise, no iterators or references are invalidated.
     *
     * ### Complexity
     *
     * At most linear in the size() of the container.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void reserve(size_type const new_cap)
    {
        data.reserve(new_cap);
    }

    /*!\brief Requests the removal of unused capacity.
     *
     * It is a non-binding request to reduce capacity() to size() and concat_capacity() to concat_size().
     * It depends on the implementation if the request is fulfilled.
     * If reallocation occurs, all iterators, including the past the end iterator, and all references to the elements
     * are invalidated. If no reallocation takes place, no iterators or references are invalidated.
     *
     * ### Complexity
     *
     * At most linear in the size() of the container.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void shrink_to_fit()
    {
        data.shrink_to_fit();
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Removes all elements from the container.
     * \returns The number of elements in the container.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    void clear() noexcept
    {
        data.clear();
    }

    /*!\brief Inserts value before position in the container.
     * \param   pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param value Element value to insert.
     * \returns     Iterator pointing to the inserted value.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ### Complexity
     *
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     */
    iterator insert(const_iterator pos, value_type const value)
    {
        return insert(pos, 1, value);
    }

    /*!\brief Inserts count copies of value before position in the container.
     * \param   pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param count Number of copies.
     * \param value Element value to insert.
     * \returns     Iterator pointing to the first element inserted, or `pos` if `count==0`.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ### Complexity
     *
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     */
    iterator insert(const_iterator pos, size_type const count, value_type const value)
    {
        auto const pos_as_num = std::distance(cbegin(), pos); // we want to insert BEFORE this position

        data.insert(data.begin() + pos_as_num, count, to_rank(value));

        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from range `[begin_it, end_it)` before position in the container.
     * \tparam begin_iterator_type Must satisfy std::ForwardIterator and
     *                             std::CommonReference<value_type_t<begin_iterator_type>, value_type>.
     * \tparam   end_iterator_type Must satisfy std::Sentinel.
     * \param[in]              pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param[in]         begin_it Begin of range to construct/assign from.
     * \param[in]           end_it End of range to construct/assign from.
     * \returns                    Iterator pointing to the first element inserted, or `pos` if `begin_it==end_it`.
     *
     * The behaviour is undefined if begin_it and end_it are iterators into `*this`.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ### Complexity
     *
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     */
    template <std::ForwardIterator begin_iterator_type, std::Sentinel<begin_iterator_type> end_iterator_type>
    iterator insert(const_iterator pos, begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires std::CommonReference<value_type_t<begin_iterator_type>, value_type>
    //!\endcond
    {
        auto const pos_as_num = std::distance(cbegin(), pos);

        auto v = ranges::iterator_range{begin_it, end_it} | seqan3::view::convert<value_type> | seqan3::view::to_rank;
        data.insert(data.begin() + pos_as_num, ranges::begin(v), ranges::end(v));

        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from initializer list before position in the container.
     * \param   pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param ilist Initializer list with values to insert.
     * \returns     Iterator pointing to the first element inserted, or `pos` if `ilist` is empty.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ### Complexity
     *
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     */
    iterator insert(const_iterator pos, std::initializer_list<value_type> const & ilist)
    {
        return insert(pos, ilist.begin(), ilist.end());
    }

    /*!\brief Removes specified elements from the container.
     * \param begin_it Begin of range to erase.
     * \param   end_it Behind the end of range to erase.
     * \returns        Iterator following the last element removed. If the iterator `pos` refers to the last element,
     *                 the end() iterator is returned.
     *
     * Invalidates iterators and references at or after the point of the erase, including the end() iterator.
     *
     * The iterator first does not need to be dereferenceable if first==end_it: erasing an empty range is a no-op.
     *
     * ### Complexity
     *
     * Linear in size().
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     */
    iterator erase(const_iterator begin_it, const_iterator end_it)
    {
        if (begin_it >= end_it) // [[unlikely]]
            return begin() + std::distance(cbegin(), end_it);

        auto const begin_it_pos = std::distance(cbegin(), begin_it);
        auto const end_it_pos = std::distance(cbegin(), end_it);

        data.erase(data.cbegin() + begin_it_pos,
                   data.cbegin() + end_it_pos);

        return begin() + begin_it_pos;
    }

    /*!\brief Removes specified elements from the container.
     * \param   pos Remove the element at pos.
     * \returns     Iterator following the last element removed. If the iterator `pos` refers to the last element,
     *              the end() iterator is returned.
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
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     */
    iterator erase(const_iterator pos)
    {
       return erase(pos, pos + 1);
    }

    /*!\brief Appends the given element value to the end of the container.
     * \param value The value to append.
     *
     * If the new size() is greater than capacity() then all iterators and references (including the past-the-end
     * iterator) are invalidated. Otherwise only the past-the-end iterator is invalidated.
     *
     * ### Complexity
     *
     * Amortised constant, worst-case linear in size().
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     */
    void push_back(value_type const value)
    {
        data.push_back(to_rank(value));
    }

    /*!\brief Removes the last element of the container.
     *
     * Calling pop_back() on an empty container is undefined. In debug mode an assertion will be thrown.
     *
     * No iterators or references except for back() and end() are invalidated.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No exception is thrown in release mode.
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void pop_back()
    {
        assert(size() > 0);
        data.pop_back();
    }

    /*!\brief Resizes the container to contain count elements.
     * \param              count The new size.
     * \throws std::length_error If count > max_size().
     * \throws    std::exception Any exception thrown by `Allocator::allocate()` (typically `std::bad_alloc`).
     *
     * Increase the size() of the vector to count.
     *
     * If the current capacity() is smaller than count, new storage is allocated and all iterators, including
     * the past-the-end iterator, and all references to the elements are invalidated.
     * Otherwise only the past-the-end iterator is invalidated.
     *
     * If the current size is greater than count, the container is reduced to its first count elements.
     * Capacity is never reduced when resizing to smaller size because that would invalidate all iterators, rather
     * than only the ones that would be invalidated by the equivalent sequence of pop_back() calls.
     *
     * ### Complexity
     *
     * At most linear in the size() of the container.
     *
     * ### Exceptions
     *
     * Only new size: Strong exception guarantee (no data is modified in case an exception is thrown).
     *
     * New default value: Basic exception guarantee, i.e. guaranteed not to leak, but container my contain bogus data
     * after exceptions is thrown.
     */
    void resize(size_type const count)
    {
        assert(count < max_size());
        data.resize(count);
    }

    /*!\copybrief resize()
     * \param value Append copies of value when resizing.
     * \copydetails resize()
     */
    void resize(size_type const count, value_type const value)
    {
        assert(count < max_size());
        data.resize(count, to_rank(value));
    }

    /*!\brief Swap contents with another instance.
     * \param rhs The other instance to swap with.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr void swap(bitcompressed_vector & rhs) noexcept
    {
        std::swap(data, rhs.data);
    }

    //!\copydoc swap()
    constexpr void swap(bitcompressed_vector && rhs) noexcept
    {
        std::swap(data, rhs.data);
    }
    //!\}

    /*!\brief Swap contents with another instance.
     * \param lhs The first instance.
     * \param rhs The other instance to swap with.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    friend constexpr void swap(bitcompressed_vector & lhs, bitcompressed_vector & rhs) noexcept
    {
        std::swap(lhs, rhs);
    }

    //!\overload
    friend constexpr void swap(bitcompressed_vector && lhs, bitcompressed_vector && rhs) noexcept
    {
        std::swap(lhs, rhs);
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(bitcompressed_vector const & rhs) const noexcept
    {
        return data == rhs.data;
    }

    constexpr bool operator!=(bitcompressed_vector const & rhs) const noexcept
    {
        return data != rhs.data;
    }

    constexpr bool operator<(bitcompressed_vector const & rhs) const noexcept
    {
        return data < rhs.data;
    }

    constexpr bool operator>(bitcompressed_vector const & rhs) const noexcept
    {
        return data > rhs.data;
    }

    constexpr bool operator<=(bitcompressed_vector const & rhs) const noexcept
    {
        return data <= rhs.data;
    }

    constexpr bool operator>=(bitcompressed_vector const & rhs) const noexcept
    {
        return data >= rhs.data;
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive_concept.
     * \param archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive_concept archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(data); //TODO: data not yet serialisable
    }
    //!\endcond
};

} // namespace seqan3
