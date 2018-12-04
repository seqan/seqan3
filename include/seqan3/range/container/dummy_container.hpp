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
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Contains the seqan3::dummy_container that throws on access but otherwise
 * behaves like a standard sequence container.
 */

#pragma once

#include <stdexcept>
#include <limits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief A dummy container that behaves like a standard sequences container except that it throws on access.
 * \tparam alphabet_type The value type of the container, must satisfy seqan3::alphabet_concept.
 * \implements seqan3::random_access_container_concept
 * \ingroup container
 *
 * This sequence is primarily used for storing a semi-alignment. We define a semi-alignment to an alignment of known
 * structure, e.g. you know it's length and the number and position of all gaps, but you do not necessarily know the
 * sequence information. The use case of a semi-alignment lies in the SAM format where alignments are stored without
 * the reference sequence information. The dummy_container let's us read in the SAM format and still reconstruct the
 * alignment with two aligned sequences (one sequence being the read/query sequence and the reference sequence being
 * represented as our dummy sequence). You can now edit the alignment or collect statistics (e.g. number of gaps)
 * with the only restriction that you **cannot access the (dummy) reference sequence.
 * See the SAM format for a detailed description of this specific use case.
 *
 * ### General Example
 *
 * The following example present some of the behaviour of the dummy_container:
 *
 * \include test/snippet/range/container/dummy_container.cpp
 *
 * The next example shows the prime use of the dummy_container in a semi-alignment.
 * Note that the dummy sequence needs to be decorated by gaps in order to fulfil
 * the seqan3::aligned_sequence_concept and become part of an alignment structure.
 *
 * \todo add example of a semi-alignment once the gap_decorators are part of seqan3.
 *
 * ### Thread safety
 *
 * This container provides no thread-safety beyond the promise given also by the STL that all
 * calls to `const` member function are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 */
template<alphabet_concept alphabet_type>
class dummy_container
{
private:
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
    //!\brief Equals alphabet_type for compatibility although actually accessing the dummy_container always throws.
    using reference         = alphabet_type;
    //!\brief Equals alphabet_type for compatibility although actually accessing the dummy_container always throws.
    using const_reference   = alphabet_type;
    //!\brief The iterator type of this container.
    //! Note that you can iterate over a dummy_container but you CANNOT access it via operator*().
    using iterator          = detail::random_access_iterator<dummy_container>;
    //!\brief The const_iterator type of this container (a random access iterator).
    //! Note that you can iterate over a dummy_container but you CANNOT access it via operator*().
    using const_iterator    = detail::random_access_iterator<dummy_container const>;
    //!\brief A signed integer type (usually std::ptrdiff_t)
    using difference_type   = ptrdiff_t;
    //!\brief An unsigned integer type (usually std::size_t)
    using size_type         = size_t;
    //!\}

    //!\cond
    // this signals to range-v3 that something is a container :|
    using allocator_type    = void;
    //!\endcond

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dummy_container() = default;
    constexpr dummy_container(dummy_container const &) = default;
    constexpr dummy_container(dummy_container &&) = default;
    constexpr dummy_container & operator=(dummy_container const &) = default;
    constexpr dummy_container & operator=(dummy_container &&) = default;
    ~dummy_container() = default;

    /*!\brief Construct from a different range.
     * \tparam    other_range_t The type of range to construct from; must satisfy std::ranges::InputRange and
     *                          std::CommonReference<value_type_t<other_range_t>, value_type>.
     * \param[in] range         The sequences to construct from.
     *
     * Note that the constructed dummy_container is of the same size as the other range but does not store any data.
     *
     * ### Complexity
     *
     * Constant if "range" is a random access range, linear in the length of "range" else.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::InputRange other_range_t>
    //!\cond
        requires has_same_value_type_v<other_range_t>
    //!\endcond
    explicit constexpr dummy_container(other_range_t && range) noexcept :
        dummy_container{std::ranges::begin(range), std::ranges::end(range)}
    {}

    /*!\brief Construct with `count` times.
     * \param[in] count Number of elements.
     * \param[in] value The initial value to be assigned (ignored).
     *
     * Note that the constructed dummy_container is of size count but does not store any data.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    dummy_container(size_type const count, value_type const value) noexcept
    {
        insert(cend(), count, value);
    }

    /*!\brief Construct from pair of iterators.
     * \tparam    begin_iterator_type Must model std::ForwardIterator and have a common reference type with value_type.
     * \tparam    end_iterator_type   Must model std::Sentinel.
     * \param[in] begin_it            Begin of range to construct from.
     * \param[in] end_it              End of range to construct from.
     *
     * Note that the constructed dummy_container is of the same size as (end_it-begin_it) but does not store any data.
     *
     * ### Complexity
     *
     * Constant if begin_it and end_it are random access iterators, linear in (end_it-begin_it) else.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ForwardIterator begin_iterator_type, std::Sentinel<begin_iterator_type> end_iterator_type>
    //!\cond
        requires std::CommonReference<value_type_t<begin_iterator_type>, value_type>
    //!\endcond
    dummy_container(begin_iterator_type begin_it, end_iterator_type end_it) noexcept
    {
        insert(cend(), begin_it, end_it);
    }

    /*!\brief Construct from `std::initializer_list`.
     * \param[in] ilist A `std::initializer_list` of value_type.
     *
     * Note that the constructed dummy_container is of the same size as the length of ilist but does not store any data.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    dummy_container(std::initializer_list<value_type> ilist) noexcept :
        dummy_container{std::begin(ilist), std::end(ilist)}
    {}

    /*!\brief Assign from `std::initializer_list`.
     * \param[in] ilist A `std::initializer_list` of value_type.
     *
     * Note that the constructed dummy_container is of the same size as the length of ilist but does not store any data.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    dummy_container & operator=(std::initializer_list<value_type> ilist) noexcept
    {
        assign(std::begin(ilist), std::end(ilist));
        return *this;
    }

    /*!\brief Assign from a different range.
     * \tparam    other_range_t The type of range to be inserted; must satisfy std::ranges::InputRange.
     * \param[in] range         The sequences to construct/assign from.
     *
     * Note that the dummy_container will have the same size as the other range but does not store any data.
     *
     * ### Complexity
     *
     * Constant if "range" is a random access range, linear in the length of "range" else.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::InputRange other_range_t>
    void assign(other_range_t && range) noexcept
    {
        dummy_container rhs{std::forward<other_range_t>(range)};
        swap(rhs);
    }

    /*!\brief Assign with `count` times `value`.
     * \param[in] count Number of elements.
     * \param[in] value The initial value to be assigned.
     *
     * Note that the dummy_container will be of length "count" but does not store any data.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    void assign(size_type const count, value_type const value) noexcept
    {
        dummy_container rhs{count, value};
        swap(rhs);
    }

    /*!\brief Assign from pair of iterators.
     * \tparam    begin_iterator_type Must model std::ForwardIterator and have a common reference type with value_type.
     * \tparam    end_iterator_type   Must satisfy std::Sentinel.
     * \param[in] begin_it            Begin of range to construct/assign from.
     * \param[in] end_it              End of range to construct/assign from.
     *
     * Note that the dummy_container will have the same size as (end_it-begin_it) but does not store any data.
     *
     * ### Complexity
     *
     * Constant if begin_it and end_it are random access iterators, linear in (end_it-begin_it) else.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ForwardIterator begin_iterator_type, std::Sentinel<begin_iterator_type> end_iterator_type>
    void assign(begin_iterator_type begin_it, end_iterator_type end_it) noexcept
    //!\cond
        requires std::CommonReference<value_type_t<begin_iterator_type>, value_type>
    //!\endcond
    {
        dummy_container rhs{begin_it, end_it};
        swap(rhs);
    }

    /*!\brief Assign from `std::initializer_list`.
     * \param[in] ilist A `std::initializer_list` of value_type.
     *
     * Note that the dummy_container will have the same size as the ilist but does not store any data.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    void assign(std::initializer_list<value_type> ilist) noexcept
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
     * Note that you can iterate over a dummy_container but cannot access the
     * value the iterator points to.
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

    /*!\name Element access (is prohibited and will throw)
     * \{
     */
    /*!\brief Throws on usage, as the dummy_container cannot not be accessed.
     * \param[in] pos Index of the element.
     * \throws std::logic_error Always, as the dummy_container cannot not be accessed.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws std::logic_error.
     */
    reference at(size_type const SEQAN3_DOXYGEN_ONLY(pos))
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
    }

    //!\copydoc at()
    const_reference at(size_type const SEQAN3_DOXYGEN_ONLY(pos)) const
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
    }

    /*!\brief Throws on usage, as the dummy_container cannot not be accessed.
     * \param[in] pos Index of the element.
     * \throws std::logic_error Always, as the dummy_container cannot not be accessed.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws std::logic_error.
     */
    reference operator[](size_type const SEQAN3_DOXYGEN_ONLY(pos))
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
    }

    //!\copydoc operator[]()
    const_reference operator[](size_type const SEQAN3_DOXYGEN_ONLY(pos)) const
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
    }

    /*!\brief Throws on usage, as the dummy_container cannot not be accessed.
     * \throws std::logic_error Always, as the dummy_container cannot not be accessed.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws std::logic_error.
     */
    reference front()
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
    }

    //!\copydoc front()
    const_reference front() const
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
    }

    /*!\brief Throws on usage, as the dummy_container cannot not be accessed.
     * \throws std::logic_error Always, as the dummy_container cannot not be accessed.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws std::logic_error.
     */
    reference back()
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
    }

    //!\copydoc back()
    const_reference back() const
    {
        throw std::logic_error{"A dummy_container cannot be accessed."};
        // end-of non-void function is never reached so the compiler allows to omit the return statement
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

    /*!\brief Returns the size of the container, i.e. std::distance(begin(), end()).
     * \returns The number of elements in the container.
     *
     * Note that the dummy container does not store any data but only the size.
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
        return size_v;
    }

    /*!\brief Returns the maximum number of elements the container is able to hold due to system or library
     * implementation limitations, i.e. std::distance(begin(), end()) for the largest container.
     * \returns The maximum possible number of elements that the container can hold.
     *
     * This value typically reflects the theoretical limit on the size of the container. At runtime, the size
     * of the container may be limited to a value smaller than max_size() by the amount of RAM available.
     *
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
        return std::numeric_limits<size_type>::max();
    }

    /*!\name Modifiers
     * \{
     */
    /*!\brief Sets the size of the container to zero.
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
        size_v = 0;
    }

    /*!\brief Increases the size of the container by one, thereby mimicking an insert without storing any data.
     * \param pos   Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param value Element value to insert.
     * \returns     Iterator pointing to the inserted value.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator insert(const_iterator pos, value_type const value) noexcept
    {
        return insert(pos, 1, value);
    }

    /*!\brief Increases the size of the container by 'count', thereby mimicking an insert without storing any data.
     * \param pos   Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param count Number of copies.
     * \param value Element value to insert (will be ignored since no data is stored).
     * \returns     Iterator pointing to the first element inserted, or `pos` if `count==0`.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator insert(const_iterator pos, size_type const count, value_type const SEQAN3_DOXYGEN_ONLY(value)) noexcept
    {
        size_v += count;
        return begin() + (pos - cbegin());
    }

    /*!\brief Increases the size of the container by (end_it - begin_it),
     *        thereby mimicking an insert without storing any data.
     * \tparam begin_iterator_type Must satisfy std::ForwardIterator and
     *                             std::CommonReference<value_type_t<begin_iterator_type>, value_type>.
     * \tparam   end_iterator_type Must satisfy std::Sentinel.
     * \param[in]              pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param[in]         begin_it Begin of range to insert from.
     * \param[in]           end_it End of range to insert from.
     * \returns                    Iterator pointing to the first element inserted, or `pos` if `begin_it==end_it`.
     *
     * ### Complexity
     *
     * Constant if begin_it and end_it are random access iterators, linear in (end_it-begin_it) else.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ForwardIterator begin_iterator_type, std::Sentinel<begin_iterator_type> end_iterator_type>
    iterator insert(const_iterator pos, begin_iterator_type begin_it, end_iterator_type end_it) noexcept
    {
        auto const diff = std::distance(begin_it, end_it);

        if (diff >= 0)
            size_v += diff;
        return begin() + (pos - cbegin());
    }

    /*!\brief Increases the size of the container by the length of the initializer list,
     *        thereby mimicking an insert without storing any data.
     * \param pos   Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param ilist Initializer list with values to insert.
     * \returns     Iterator pointing to the first element inserted, or `pos` if `ilist` is empty.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator insert(const_iterator pos, std::initializer_list<value_type> const & ilist) noexcept
    {
        return insert(pos, ilist.begin(), ilist.end());
    }

    /*!\brief Decreases the size of the container by (end_it - begin_it),
     *        thereby mimicking an erase without storing any data.
     * \param begin_it Begin of range to erase.
     * \param end_it   Behind the end of range to erase.
     * \returns        Iterator following the last element removed. If the iterator `pos` refers to the last element,
     *                 the end() iterator is returned.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator erase(const_iterator begin_it, const_iterator end_it) noexcept
    {
        if (begin_it >= end_it) // [[unlikely]]
            return begin() + (end_it - cbegin());

        size_v -= (end_it - begin_it);

        return begin() + (begin_it - cbegin());
    }

    /*!\brief Decreases the size of the container by one, thereby mimicking an erase without storing any data.
     * \param   pos Remove the element at pos.
     * \returns     Iterator following the last element removed. If the iterator `pos` refers to the last element,
     *              the end() iterator is returned.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator erase(const_iterator pos) noexcept
    {
       return erase(pos, pos + 1);
    }

    /*!\brief Increases the size of the container by one, thereby mimicking a puch_back without storing any data.
     * \param value The value to append (ignored).
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    void push_back(value_type const value) noexcept
    {
        ++size_v;
    }

    /*!\brief Decreases the size of the container by one, thereby mimicking a pop_back without storing any data.
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
     * No-throw guarantee.
     */
    void pop_back() noexcept
    {
        assert(size_v > 0);
        --size_v;
    }

    /*!\brief Resizes the container to the size of 'count' elements.
     * \param              count The new size.
     * \throws std::length_error If count > max_size().
     * \throws    std::exception Any exception thrown by `Allocator::allocate()` (typically `std::bad_alloc`).
     *
     * Set the size() of the vector to count.
     *
     * Since the dummy_container does not store any data, no memory allocation or deallocation is performed.
     * If count < size(), all iterators after count are invalidated.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    void resize(size_type const count) noexcept
    {
        assert(count < max_size());
        size_v = count;
    }

    /*!\copybrief resize()
     * \param value Append copies of value when resizing (ignored).
     * \copydetails resize()
     */
    void resize(size_type const count, value_type const SEQAN3_DOXYGEN_ONLY(value))
    {
        resize(count);
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
    constexpr void swap(dummy_container & rhs) noexcept
    {
        std::swap(size_v, rhs.size_v);
    }

    //!\copydoc swap()
    constexpr void swap(dummy_container && rhs) noexcept
    {
        std::swap(size_v, rhs.size_v);
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
    friend constexpr void swap(dummy_container & lhs, dummy_container & rhs) noexcept
    {
        std::swap(lhs, rhs);
    }

    //!\overload
    friend constexpr void swap(dummy_container && lhs, dummy_container && rhs) noexcept
    {
        std::swap(lhs, rhs);
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(dummy_container const & rhs) const noexcept
    {
        return size_v == rhs.size_v;
    }

    constexpr bool operator!=(dummy_container const & rhs) const noexcept
    {
        return size_v != rhs.size_v;
    }

    constexpr bool operator<(dummy_container const & rhs) const noexcept
    {
        return size_v < rhs.size_v;
    }

    constexpr bool operator>(dummy_container const & rhs) const noexcept
    {
        return size_v > rhs.size_v;
    }

    constexpr bool operator<=(dummy_container const & rhs) const noexcept
    {
        return size_v <= rhs.size_v;
    }

    constexpr bool operator>=(dummy_container const & rhs) const noexcept
    {
        return size_v >= rhs.size_v;
    }
    //!\}

private:
    size_type size_v{0};
};


} // namespace seqan3
