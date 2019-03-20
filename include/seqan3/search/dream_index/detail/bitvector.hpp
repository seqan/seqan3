// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::detail::bitvector.
 */

#pragma once

#include <sdsl/bit_vectors.hpp>

// This file will take care of
//  * Uncompressed bitvectors
//  * Compressed bitvectors
//  * Chunked bitvector (chunks can be (un)-compressed)

// Basically a wrapper that provides a unified interface for all specialisations
// => How to handle chunked?

// This should not be in detail since it is used for traits.
namespace seqan3
{

    //!\brief Tag for the uncompressed vector.
    struct uncompressed {};
    //!\brief Tag for the compressed vector.
    struct compressed {};

    /*!\brief Tag for the chunked vector.
     * \tparam chunks The number of chunks.
     */
    template <size_t chunks>
    struct chunked
    {
        //!\brief The number of chunks.
        static constexpr size_t value{chunks};
    };

} // namespace seqan3

namespace seqan3::detail
{

template <typename strategy>
class bitvector {};

//!\brief The uncompressed bitvector.
template <>
class bitvector<uncompressed> // Implements reservable_container_concept
{
private:
    //!\brief The underlying sdsl vector implementation.
    using sdsl_vector_type = sdsl::bit_vector;
    //!\brief The underlying data structure.
    sdsl_vector_type data_value;

public:
    /*!\name Member types
     * \{
     */
    //!\brief An unsigned integer type (usually std::size_t).
    //!\hideinitializer
    using value_type = sdsl_vector_type::value_type;

    //!\brief The reference type of this container.
    //!\hideinitializer
    using reference = sdsl_vector_type::reference;

    //!\brief The const reference type of this container.
    //!\hideinitializer
    using const_reference = sdsl_vector_type::const_reference;

    //!\brief An unsigned integer type (usually std::size_t).
    //!\hideinitializer
    using size_type = sdsl_vector_type::size_type;

    //!\brief The iterator type of this container.
    //!\hideinitializer
    using iterator = sdsl_vector_type::iterator;

    //!\brief The const reference type of this container.
    //!\hideinitializer
    using const_iterator = sdsl_vector_type::const_iterator;
    //!\}

    //!\cond
    // this signals to range-v3 that something is a container :|
    using allocator_type = void;
    //!\endcond

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bitvector() = default;                              //!< Default constructor.
    bitvector(bitvector const &) = default;             //!< Copy constructor.
    bitvector & operator=(bitvector const &) = default; //!< Move constructor.
    bitvector(bitvector &&) = default;                  //!< Copy assignment.
    bitvector & operator=(bitvector &&) = default;      //!< Move assignment.
    ~bitvector() = default;                             //!< Destructor.

    /*!\brief Construct/assign with `count` times `value`.
     * \param count Number of elements.
     * \param value The initial value to be assigned.
     *
     * \par Complexity
     *
     * In \f$O(count*value)\f$.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    bitvector(size_type const count, value_type value = 0)
    {
        data_value = sdsl::bit_vector(count, value);
    }

    /*!\brief Construct/assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy std::ForwardIterator.
     * \tparam end_iterator_type Must satisfy std::SizedSentinel.
     * \param begin_it begin of range to construct/assign from.
     * \param end_it   end of range to construct/assign from.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of the ranges between `begin_it` and `end_it`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ForwardIterator begin_iterator_type, std::SizedSentinel<begin_iterator_type> end_iterator_type>
    bitvector(begin_iterator_type begin_it, end_iterator_type end_it)
    {
        insert(cend(), begin_it, end_it);
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \param ilist an `std::initializer_list` of `value_type`.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of the ranges in `ilist`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    bitvector(std::initializer_list<value_type> const & ilist)
    {
        data_value.assign(ilist);
    }
    //!\}

    /*!\brief Return the i-th element.
     * \param i The element to retrieve.
     * \returns A reference to the i-th element.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (never modifies data)..
     */
    reference operator[](size_type const i) noexcept
    {
        return data_value[i]; // sdsl asserts in debug mode
    }

    //!\copydoc operator[]()
    const_reference operator[](size_type const i) const noexcept
    {
        return data_value[i]; // sdsl asserts in debug mode
    }

    /*!\brief Return the i-th element.
     * \param i The element to retrieve.
     * \throws std::out_of_range If you access an element behind the last.
     * \returns A reference to the i-th element.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (never modifies data)..
     */
    reference at(size_type const i)
    {
        if (i >= size()) // [[unlikely]]
            throw std::out_of_range{"Trying to access element behind the last in bitvector."};
        return (*this)[i];
    }

    //!\copydoc at()
    const_reference at(size_type const i) const
    {
        if (i >= size()) // [[unlikely]]
            throw std::out_of_range{"Trying to access element behind the last in bitvector."};
        return (*this)[i];
    }

    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return iterator{&data_value};
    }

    //!\copydoc begin()
    const_iterator begin() const noexcept
    {
        return const_iterator{&data_value};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
    {
        return const_iterator{&data_value};
    }

    /*!\brief Returns an iterator to the element following the last element of the container.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    iterator end() noexcept
    {
        return iterator{&data_value, size()};
    }

    //!\copydoc end()
    const_iterator end() const noexcept
    {
        return const_iterator{&data_value, size()};
    }

    //!\copydoc end()
    const_iterator cend() const noexcept
    {
        return const_iterator{&data_value, size()};
    }

    //!\name Comparison operators
    //!\{
    /*!\brief Test for equality.
     * \param rhs Bitvector to compare to.
     * \returns `true` if equal, `false` otherwise.
     */
    bool operator==(bitvector const & rhs) const noexcept
    {
        return data() == rhs.data();
    }

    /*!\brief Test for inequality.
     * \param rhs Bitvector to compare to.
     * \returns `true` if unequal, `false` otherwise.
     */
    bool operator!=(bitvector const & rhs) const noexcept
    {
        return data() != rhs.data();
    }

    /*!\brief Test for less than.
     * \param rhs Bitvector to compare to.
     * \returns `true` if this is smaller than rhs, `false` otherwise.
     */
    bool operator<(bitvector const & rhs) const noexcept
    {
        return data() < rhs.data();
    }

    /*!\brief Test for greater than.
     * \param rhs Bitvector to compare to.
     * \returns `true` if this is greater than rhs, `false` otherwise.
     */
    bool operator>(bitvector const & rhs) const noexcept
    {
        return data() > rhs.data();
    }

    /*!\brief Test for less than or equal.
     * \param rhs Bitvector to compare to.
     * \returns `true` if this is less than or equal to rhs, `false` otherwise.
     */
    bool operator<=(bitvector const & rhs) const noexcept
    {
        return data() <= rhs.data();
    }

    /*!\brief Test for greater than or equal.
     * \param rhs Bitvector to compare to.
     * \returns `true` if this is greater than or equal to rhs, `false` otherwise.
     */
    bool operator>=(bitvector const & rhs) const noexcept
    {
        return data() >= rhs.data();
    }
    //!\}

    /*!\brief Swap contents with another instance.
     * \param rhs The other instance to swap with.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    void swap(bitvector & rhs) noexcept
    {
        std::swap(data_value, rhs.data_value);
    }

    //!\copydoc swap()
    void swap(bitvector && rhs) noexcept
    {
        std::swap(data_value, rhs.data_value);
    }

    /*!\brief Assign a value to a position in the bitvector.
     * \param index The position to assign to.
     * \param value The value to assign.
     *
     * \attention This method will set as many bits as there are in `size_type`, i.e. usually 64 bits.
     *
     * \par Complexity
     *
     * Constant.
     */
    void set_int(size_type const index, size_type const value)
    {
        data_value.set_int(index, value);
    }

    /*!\brief Returns the bits at a position.
     * \param index The position to read.
     * \returns A size_t where each bit represent a bit starting at index of the bitvector.
     *
     * \attention This method will read as many bits as there are in `size_type`, i.e. usually 64 bits.
     *
     * \par Complexity
     *
     * Constant.
     */
    size_t get_int(size_type const index) const noexcept
    {
        return data_value.get_int(index);
    }

    /*!\brief Returns the number of elements in the container, i.e. std::distance(begin(), end()).
     * \returns The number of elements in the container.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
    {
        return data_value.size();
    }

    /*!\brief Returns the maximum number of elements the container is able to hold due to system or library
     * implementation limitations, i.e. std::distance(begin(), end()) for the largest container.
     * \returns The number of elements in the container.
     *
     * This value typically reflects the theoretical limit on the size of the container. At runtime, the size
     * of the container may be limited to a value smaller than max_size() by the amount of RAM available.
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    size_type max_size() const noexcept
    {
        return data_value.max_size();
    }

    /*!\brief Checks whether the container is empty.
     * \returns `true` if the container is empty, `false` otherwise.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    bool empty() const noexcept
    {
        return data_value.empty();
    }

    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns sdsl::bit_vector.
     *
     * This exact representation of the data is implementation defined. Do not rely on it for API stability.
     */
    sdsl_vector_type & data()
    {
        return data_value;
    }

    //!\copydoc data()
    sdsl_vector_type const & data() const
    {
        return std::as_const(data_value);
    }

    /*!\brief Return the first element. Calling front on an empty container is undefined.
     * \returns A reference to the first element.
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (never modifies data).
     */
    reference front() noexcept
    {
        assert(size() > 0);
        return data_value.front();
    }

    //!\copydoc front()
    const_reference front() const noexcept
    {
        assert(size() > 0);
        return data_value.front();
    }

    /*!\brief Return the last element.
     * \returns A reference to the last element.
     *
     * Calling back on an empty container is undefined. In debug mode an assertion checks the size of the container.
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (never modifies data)..
     */
    reference back() noexcept
    {
        assert(size() > 0);
        return data_value.back();
    }

    //!\copydoc back()
    const_reference back() const noexcept
    {
        assert(size() > 0);
        return data_value.back();
    }

    /*!\brief Inserts value before position in the container.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param value Element value to insert.
     * \returns Iterator pointing to the inserted value.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \par Complexity
     *
     * Constant plus linear in the distance between pos and end of the container.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    iterator insert(const_iterator pos, value_type const value)
    {
        return data_value.insert(pos, value);
    }

    /*!\brief Inserts value before position in the container.
     * \param count Number of copies.
     * \returns Iterator pointing to the first element inserted, or pos if `count==0`.
     * \copydetails insert(const_iterator pos, value_type const value)
     */
    iterator insert(const_iterator pos, size_type const count, value_type const value)
    {
        return data_value.insert(pos, count, value);
    }

    /*!\brief Inserts elements from range `[first, last)` before position in the container.
     * \tparam begin_iterator_type Must satisfy std::ForwardIterator.
     * \tparam end_iterator_type Must satisfy std::SizedSentinel.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param first Begin of range to insert.
     * \param last Behind the end of range to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `first==last`.
     *
     * The behaviour is undefined if first and last are iterators into `*this`.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \par Complexity
     *
     * Constant plus linear in the distance between pos and end of the container.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     *
     */
    template <std::ForwardIterator begin_iterator_type, typename end_iterator_type>
    iterator insert(const_iterator pos, begin_iterator_type first, end_iterator_type last)
    //!\cond
        requires std::SizedSentinel<end_iterator_type, begin_iterator_type>
    //!\endcond
    {
        return data_value.insert(pos, first, last);
    }

    /*!\brief Inserts elements from initializer list before position in the container.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param ilist Initializer list with values to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `ilist` is empty.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \par Complexity
     *
     * Constant plus linear in the distance between pos and end of the container.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    iterator insert(const_iterator pos, std::initializer_list<value_type> const & ilist)
    {
        return data_value.insert(pos, ilist);
    }

    /*!\brief Construct/assign with `count` times `value`.
     * \param count Number of elements.
     * \param value The initial value to be assigned.
     *
     * \par Complexity
     *
     * In \f$O(count*value)\f$.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void assign(size_type const count, value_type const value)
    {
        data_value.assign(count, value);
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \param ilist an `std::initializer_list` of `value_type`.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of the elements in `ilist`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void assign(std::initializer_list<value_type> const & ilist)
    {
        data_value.assign(ilist);
    }

    /*!\brief Construct/assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy std::ForwardIterator.
     * \tparam end_iterator_type Must satisfy std::SizedSentinel.
     * \param begin_it begin of range to construct/assign from.
     * \param end_it end of range to construct/assign from.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of the ranges between `begin_it` and `end_it`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ForwardIterator begin_iterator_type, typename end_iterator_type>
    void assign(begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires std::SizedSentinel<end_iterator_type, begin_iterator_type>
    //!\endcond
    {
        data_value.assign(begin_it, end_it);
    }

    /*!\brief Removes specified elements from the container.
     * \param first Begin of range to erase.
     * \param last Behind the end of range to erase.
     * \returns Iterator pointing to the first element inserted, or pos if `first==last`.
     *
     * Invalidates iterators and references at or after the point of the erase, including the end() iterator.
     *
     * The iterator first does not need to be dereferenceable if first==last: erasing an empty range is a no-op.
     *
     * \par Complexity
     *
     * Linear in concat_size().
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    iterator erase(const_iterator first, const_iterator last)
    {
        return data_value.erase(first, last);
    }

    /*!\brief Removes specified elements from the container.
     * \param pos Remove the element at pos.
     * \returns Iterator pointing to the first element inserted, or pos if `first==last`.
     *
     * Invalidates iterators and references at or after the point of the erase, including the end() iterator.
     *
     * The iterator pos must be valid and dereferenceable. Thus the end() iterator (which is valid, but is not
     * dereferencable) cannot be used as a value for pos.
     *
     * \par Complexity
     *
     * Linear in concat_size().
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    iterator erase(const_iterator pos)
    {
        return data_value.erase(pos);
    }

    /*!\brief Appends the given element value to the end of the container.
     * \tparam rng_type The type of range to be inserted.
     * \param value The value to append.
     *
     * If the new size() is greater than capacity() then all iterators and references (including the past-the-end
     * iterator) are invalidated. Otherwise only the past-the-end iterator is invalidated.
     *
     * \par Complexity
     *
     * Amortised constant.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    void push_back(value_type const value)
    {
        data_value.push_back(value);
    }

    /*!\brief Removes the last element of the container.
     *
     * Calling pop_back on an empty container is undefined. In debug mode an assertion will be thrown.
     *
     * No iterators or references except for back() and end() are invalidated.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No exception is thrown in release mode.
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void pop_back()
    {
        data_value.pop_back();
    }

    /*!\brief Removes all elements from the container.
    *
    * \par Complexity
    *
    * Constant.
    *
    * \par Exceptions
    *
    * No-throw guarantee.
     */
    void clear()
    {
        data_value.clear();
    }

    /*!\brief Requests the removal of unused capacity.
     *
     * It is a non-binding request to reduce capacity() to size().
     * It depends on the implementation if the request is fulfilled.
     * If reallocation occurs, all iterators, including the past the end iterator, and all references to the elements
     * are invalidated. If no reallocation takes place, no iterators or references are invalidated.
     *
     * \par Complexity
     *
     * At most linear in the size() of the container.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void shrink_to_fit()
    {
        data_value.shrink_to_fit();
    }

    /*!\brief Increase the capacity to a value that's greater or equal to new_cap.
     * \param new_cap The new capacity.
     * \throws std::bad_alloc If memory allocation failed.
     *
     * Increase the capacity of the vector to a value that's greater or equal to new_cap.
     * If new_cap is greater than the current capacity(), new storage is allocated, otherwise the method does nothing.
     * If new_cap is greater than capacity(), all iterators, including the past-the-end iterator, and all references
     * to the elements are invalidated. Otherwise, no iterators or references are invalidated.
     *
     * \par Complexity
     *
     * At most linear in the size() of the container.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void reserve(size_type const new_cap)
    {
        data_value.reserve(new_cap);
    }

    /*!\brief Resizes the container to contain count elements.
     * \param count The new size.
     * \throws std::bad_alloc If memory allocation failed.
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
     * \par Complexity
     *
     * At most linear in the size() of the container.
     *
     * \par Exceptions
     *
     * Only new size: Strong exception guarantee (no data is modified in case an exception is thrown). [only new size]
     *
     * New default value: Basic exception guarantee, i.e. guaranteed not to leak, but container my contain bogus data
     * after exceptions is thrown.
     */
    void resize(size_type const count)
    {
        data_value.resize(count);
    }

    /*!\copybrief resize()
     * \tparam rng_type The type of range to be inserted.
     * \param value Instead of appending empty containers, append copies of value.
     * \copydetails resize()
     */
    void resize(size_type const count, value_type const value)
    {
        data_value.resize(count, value);
    }

    /*!\brief Returns the number of elements that the container has currently allocated space for.
     * \returns The capacity of the currently allocated storage.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    size_type capacity() const noexcept
    {
        return data_value.capacity();
    }
};

} // namespace seqan3::detail

namespace seqan3
{

    //!\overload
    template <typename strategy>
    void swap(detail::bitvector<strategy> lhs, detail::bitvector<strategy> rhs)
    {
        lhs.swap(rhs);
    }

} // namespace seqan3
