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

/*!\file
 * \brief Provides seqan3::concatenated_sequences.
 * \ingroup container
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>
#include <vector>

#include <range/v3/view/bounded.hpp>
#include <range/v3/view/const.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <range/v3/view/slice.hpp>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/metafunction.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/types/vector.hpp>
#endif

namespace seqan3
{

/*!\brief Container that stores sequences concatenated internally.
 * \tparam inner_type The type of sequences that will be stored. Must satisfy seqan3::reservable_sequence_concept.
 * \tparam data_delimiters_type A container that stores the begin/end positions in the inner_type. Must be
 * seqan3::reservable_sequence_concept and have inner_type's size_type as value_type.
 * \implements seqan3::reservable_sequence_concept
 * \ingroup container
 *
 * This class may be used whenever you would usually use `std::vector<std::vector<some_alphabet>>` or
 * `std::vector<std::string>`, i.e. whenever you have a collection of sequences. It is the spiritual successor of
 * the `StringSet<TString, Owner<ConcatDirect>>` from SeqAn2.
 *
 * It saves all of the member sequences inside one concatenated sequence internally. If you access an element,
 * you instead get a view on the internal string as a proxy. This has the following
 * advantages:
 *
 * * Better cache locality when parsing the sequences linearly (and often also on random access).
 * * Constant time access to the concatenation of the sequences via concat().
 * * This access is also writable so that certain transformations can be done globally, instead of element-wise.
 * * Also direct access to the delimiters via data() [this is used by some algorithms].
 *
 * The disadvantages are:
 *
 * * Slower inserts and erases because the entire concatenation might have to be copied.
 * * No emplace operations.
 * * Modifying elements is limited to operations on elements of that element, i.e. you can change a character,
 * but you can't assign a new member sequence to an existing position.
 *
 * \par Example
 *
 * ```cpp
 * concatenated_sequences<dna4_vector> concat1{"ACGT"_dna4, "GAGGA"_dna4};
 * std::cout << concat1[0] << '\n'; // "[A,C,G,T]"
 *
 * std::vector<dna4_vector> concat2{"ACTA"_dna4, "AGGA"_dna4};
 *
 * concat1 = concat2;               // you can assign from other ranges
 *
 * concat2[0] = "ATTA"_dna4;        // this works for vector of vector
 * //concat1[0] = "ATTA"_dna4;      // but not on concatenated_sequences
 *
 * concat1[0][1] = dna4::T;         // this, however, does
 * std::cout << concat1[0] << '\n'; // "[A,T,T,A]"
 *
 *
 * // if you know that you will be adding a thousand vectors of length thousand:
 * concat1.reserve(1'000);
 * concat1.concat_reserve(1'000 * 1'000);
 * while (...)
 * {
 *     // ...
 *     concat1.push_back(vector_of_length1000);
 * }
 * ```
 *
 * \par Exceptions
 *
 * Whenever a strong exception guarantee is given for this class, it presumes that
 * `std::is_nothrow_move_constructible<typename inner_type::value_type>` otherwise only basic exception safety can
 * be assumed.
 *
 * \par Thread safety
 *
 * This container provides no thread-safety beyond the promise given also by the STL that all
 * calls to `const` member function are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 */
template <typename inner_type,
          typename data_delimiters_type = std::vector<typename inner_type::size_type>>
//!\cond
    requires reservable_sequence_concept<std::remove_reference_t<inner_type>> &&
             reservable_sequence_concept<std::remove_reference_t<data_delimiters_type>> &&
             std::is_same_v<ranges::size_type_t<inner_type>, ranges::value_type_t<data_delimiters_type>>
//!\endcond
class concatenated_sequences
{
    static_assert(!detail::sequence_concept_modified_by_const_iterator_bug<inner_type>,
                  "KNOWN BUG: inner_type = std::basic_string<> is not working "
                  "for the ubuntu::ppa version of gcc7, because of a faulty STL version. ");
protected:
    //!\privatesection
    //!\brief Where the concatenation is stored.
    std::decay_t<inner_type> data_values;
    //!\brief Where the delimiters are stored; begins with 0, has size of size() + 1.
    data_delimiters_type data_delimiters{0};

public:
    //!\publicsection
    /*!\name Member types
     * \{
     */
    //!\brief == inner_type.
    //!\hideinitializer
    using value_type = std::decay_t<inner_type>;

    //!\brief A proxy of type ranges::view::slice that represents the range on the concatenated vector.
    //!\hideinitializer
    using reference = decltype(data_values | ranges::view::slice(0, 1));

    //!\brief An immutable proxy of type ranges::view::slice that represents the range on the concatenated vector.
    //!\hideinitializer
    using const_reference = decltype(std::as_const(data_values) | ranges::view::slice(0, 1) | ranges::view::const_);

    //!\brief The iterator type of this container (a random access iterator).
    //!\hideinitializer
    using iterator = detail::random_access_iterator<concatenated_sequences>;

    //!\brief The const iterator type of this container (a random access iterator).
    //!\hideinitializer
    using const_iterator = detail::random_access_iterator<concatenated_sequences const>;

    //!\brief A signed integer type (usually std::ptrdiff_t)
    //!\hideinitializer
    using difference_type = ranges::difference_type_t<data_delimiters_type>;

    //!\brief An unsigned integer type (usually std::size_t)
    //!\hideinitializer
    using size_type = ranges::size_type_t<data_delimiters_type>;
    //!\}

    /* rule of six */
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructors.
    concatenated_sequences() = default;
    //!\brief Default constructors.
    constexpr concatenated_sequences(concatenated_sequences const &) = default;
    //!\brief Default constructors.
    constexpr concatenated_sequences(concatenated_sequences &&) = default;
    //!\brief Default constructors.
    constexpr concatenated_sequences & operator=(concatenated_sequences const &) = default;
    //!\brief Default constructors.
    constexpr concatenated_sequences & operator=(concatenated_sequences &&) = default;
    //!\brief Default constructors.
    ~concatenated_sequences() = default;

    /*!\brief Construct/assign from a different range.
     * \tparam rng_of_rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with
     * `concatenated_sequences`.
     * \param rng_of_rng The sequences to construct/assign from.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of `rng_of_rng`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <typename rng_of_rng_type>
    concatenated_sequences(rng_of_rng_type && rng_of_rng)
    //!\cond
        requires input_range_concept<std::decay_t<rng_of_rng_type>> &&
                 compatible_concept<rng_of_rng_type, concatenated_sequences>
    //!\endcond
    {
        if constexpr (sized_range_concept<std::decay_t<rng_of_rng_type>>)
            data_delimiters.reserve(ranges::size(rng_of_rng) + 1);

        for (auto && val : rng_of_rng)
        {
            data_values.insert(data_values.end(), val.begin(), val.end());
            data_delimiters.push_back(data_delimiters.back() + val.size());
        }
    }

    /*!\brief Construct/assign with `count` times `value`.
     * \tparam rng_type The type of range to be inserted; `rng_type` and `value_type` must satisfy
     * seqan3::compatible_concept.
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
    template <typename rng_type>
    concatenated_sequences(size_type const count, rng_type && value)
    //!\cond
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    //!\endcond
    {
        // TODO SEQAN_UNLIKELY
        if (count == 0)
            return;

        insert(cend(), count, std::forward<rng_type>(value));
    }

    /*!\brief Construct/assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy seqan3::forward_iterator_concept and satisfy seqan3::compatible_concept
     * with `value_type`.
     * \tparam end_iterator_type Must satisfy seqan3::sized_sentinel_concept.
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
    template <typename begin_iterator_type, typename end_iterator_type>
    concatenated_sequences(begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires forward_iterator_concept<begin_iterator_type> &&
                 compatible_concept<begin_iterator_type, concatenated_sequences>
                 //&& sized_sentinel_concept<end_iterator_type, begin_iterator_type>
    //!\endcond
    {
        insert(cend(), begin_it, end_it);
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \tparam rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with `value_type`.
     * \param ilist an `std::initializer_list` of `rng_type`.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of the ranges in `ilist`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <typename rng_type = value_type>
    concatenated_sequences(std::initializer_list<rng_type> ilist)
    //!\cond
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    //!\endcond
    {
        assign(std::begin(ilist), std::end(ilist));
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \tparam rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with `value_type`.
     * \param ilist an `std::initializer_list` of `rng_type`.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of the elements in `ilist`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <typename rng_type>
    concatenated_sequences & operator=(std::initializer_list<rng_type> ilist)
    //!\cond
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    //!\endcond
    {
        assign(std::begin(ilist), std::end(ilist));
        return *this;
    }

    /*!\brief Construct/assign from a different range.
     * \tparam rng_of_rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with
     * `concatenated_sequences`.
     * \param rng_of_rng The sequences to construct/assign from.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of `rng_of_rng`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <typename rng_of_rng_type>
    void assign(rng_of_rng_type && rng_of_rng)
    //!\cond
        requires input_range_concept<std::decay_t<rng_of_rng_type>> &&
                 compatible_concept<rng_of_rng_type, concatenated_sequences>
    //!\endcond
    {
        concatenated_sequences rhs{std::forward<rng_of_rng_type>(rng_of_rng)};
        swap(rhs);
    }

    /*!\brief Construct/assign with `count` times `value`.
     * \tparam rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with `value_type`.
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
    template <typename rng_type>
    void assign(size_type const count, rng_type && value)
    //!\cond
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    //!\endcond
    {
        concatenated_sequences rhs{count, value};
        swap(rhs);
    }

    /*!\brief Construct/assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy seqan3::forward_iterator_concept and satisfy seqan3::compatible_concept
     * with `value_type`.
     * \tparam end_iterator_type Must satisfy seqan3::sized_sentinel_concept.
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
    template <typename begin_iterator_type, typename end_iterator_type>
    void assign(begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires forward_iterator_concept<begin_iterator_type> &&
                 compatible_concept<begin_iterator_type, concatenated_sequences>
                 //&& sized_sentinel_concept<end_iterator_type, begin_iterator_type>
    //!\endcond
    {
        concatenated_sequences rhs{begin_it, end_it};
        swap(rhs);
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \tparam rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with `value_type`.
     * \param ilist an `std::initializer_list` of `rng_type`.
     *
     * \par Complexity
     *
     * Linear in the cumulative size of the elements in `ilist`.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <typename rng_type = value_type>
    void assign(std::initializer_list<rng_type> ilist)
    //!\cond
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    //!\endcond
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
    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \throws std::out_of_range If you access an element behind the last.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
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
        //TODO add SEQAN_UNLIKELY
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in concatenated_sequences."};
        return (*this)[i];
    }

    //!\copydoc at()
    const_reference at(size_type const i) const
    {
        //TODO add SEQAN_UNLIKELY
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in concatenated_sequences."};
        return (*this)[i];
    }

    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
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
    reference operator[](size_type const i)
    {
        assert(i < size());
        return data_values | ranges::view::slice(data_delimiters[i], data_delimiters[i+1]);
    }

    //!\copydoc operator[]()
    const_reference operator[](size_type const i) const
    {
        assert(i < size());
        return data_values | ranges::view::slice(data_delimiters[i], data_delimiters[i+1])
                           | ranges::view::const_;
    }

    /*!\brief Return the first element as a view. Calling front on an empty container is undefined.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
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
    reference front()
    {
        assert(size() > 0);
        return (*this)[0];
    }

    //!\copydoc front()
    const_reference front() const
    {
        assert(size() > 0);
        return (*this)[0];
    }

    /*!\brief Return the last element as a view.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
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
    reference back()
    {
        assert(size() > 0);
        return (*this)[size()-1];
    }

    //!\copydoc back()
    const_reference back() const
    {
        assert(size() > 0);
        return (*this)[size()-1];
    }

    /*!\brief Return the concatenation of all members.
     * \returns A ranges::view proxy on the concatenation of underlying sequences.
     *
     * This is a safe way of accessing the internal concatenated representation, i.e. you cannot do operations
     * that would invalidate this container (like insert or resize), but you can write to the individual positions.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (never modifies data).
     */
    reference concat()
    {
        return data_values | ranges::view::slice(static_cast<size_type>(0), concat_size());
    }

    //!\copydoc concat()
    const_reference concat() const
    {
        return data_values | ranges::view::slice(static_cast<size_type>(0), concat_size()) | ranges::view::const_;
    }

    /*!\brief Provides direct, unsafe access to underlying data structures.
     * \returns An std::pair of the concatenated sequences and the delimiter string.
     *
     * This exact representation of the data is implementation defined. Do not rely on it for API stability.
     */
    std::pair<decltype(data_values) &, decltype(data_delimiters) &> data()
    {
        return {data_values, data_delimiters};
    }

    //!\copydoc data()
    std::pair<decltype(data_values) const &, decltype(data_delimiters) const &> data() const
    {
        return {std::as_const(data_values), std::as_const(data_delimiters)};
    }
    //!\}

    /*!\name Capacity
     * \{
     */
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
        return size() == 0;
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
        return data_delimiters.size() - 1;
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
        return data_delimiters.max_size() - 1;
    }

    /*!\brief Returns the number of elements that the container has currently allocated space for.
     * \returns The capacity of the currently allocated storage.
     *
     * \attention
     *
     * This does not operate on underlying concat container, see concat_capacity().
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
        return data_delimiters.capacity();
    }

    /*!\brief Increase the capacity to a value that's greater or equal to new_cap.
     * \param new_cap The new capacity.
     * \throws std::length_error If new_cap > max_size().
     * \throws std::exception Any exception thrown by `Allocator::allocate()` (typically `std::bad_alloc`).
     *
     * Increase the capacity of the vector to a value that's greater or equal to new_cap.
     * If new_cap is greater than the current capacity(), new storage is allocated, otherwise the method does nothing.
     * If new_cap is greater than capacity(), all iterators, including the past-the-end iterator, and all references
     * to the elements are invalidated. Otherwise, no iterators or references are invalidated.
     *
     * \attention
     *
     * This does not operate on underlying concat container, see concat_reserve().
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
        data_delimiters.reserve(new_cap + 1);
    }

    /*!\brief Requests the removal of unused capacity.
     *
     * It is a non-binding request to reduce capacity() to size() and concat_capacity() to concat_size().
     * It depends on the implementation if the request is fulfilled.
     * If reallocation occurs, all iterators, including the past the end iterator, and all references to the elements
     * are invalidated. If no reallocation takes place, no iterators or references are invalidated.
     *
     * \attention
     *
     * This effects both underlying data structures.
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
        data_values.shrink_to_fit();
        data_delimiters.shrink_to_fit();
    }
    //!\}

    /*!\name Capacity (concat)
     * \{
     */
    /*!\brief Returns the cumulative size of all elements in the container.
     * \returns The cumulative size of elements in the container.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    size_type concat_size() const noexcept
    {
        return data_values.size();
    }

    /*!\brief Returns the concatenated size the container has currently allocated space for.
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
    size_type concat_capacity() const noexcept
    {
        return data_values.capacity();
    }

    /*!\brief Increase the concat_capacity() to a value that's greater or equal to new_cap.
     * \param new_cap The new capacity.
     * \throws std::length_error If new_cap > max_size().
     * \throws std::exception Any exception thrown by `Allocator::allocate()` (typically `std::bad_alloc`).
     *
     * Increase the capacity of the underlying concatenated sequence to a value that's greater or equal to new_cap.
     * If new_cap is greater than the current concat_capacity(), new storage is allocated, otherwise the method does
     * nothing. If new_cap is greater than concat_capacity(), all iterators, including the past-the-end iterator, and
     * all references to the elements are invalidated. Otherwise, no iterators or references are invalidated.
     *
     * \par Complexity
     *
     * At most linear in the concat_size() of the container.
     *
     * \par Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    void concat_reserve(size_type const new_cap)
    {
        data_values.reserve(new_cap);
    }
    //!\}


    /*!\name Modifiers
     * \{
     */
    /*!\brief Removes all elements from the container.
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
    void clear() noexcept
    {
        data_values.clear();
        data_delimiters.clear();
        data_delimiters.push_back(0);
    }

    /*!\brief Inserts value before position in the container.
     * \tparam rng_type The type of range to be inserted; must satisfy seqan3::forward_range_concept
     * and have the same `value_type` as `value_type` (i.e. `value_type`'s `value_type`!).
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
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     *
     * \par Example
     *
     * ```cpp
     * concatenated_sequences<dna4_vector> foobar;
     * foobar.insert(foobar.end(), "ACGT"_dna4);
     * std::cout << foobar[0] << '\n'; // [A, C, G, T]
     * ```
     */
    template <typename rng_type>
    iterator insert(const_iterator pos, rng_type && value)
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    {
        return insert(pos, 1, std::forward<rng_type>(value));
    }
    // no specialisation for temporaries, since we have to copy anyway

    /*!\brief Inserts count copies of value before position in the container.
     * \tparam rng_type The type of range to be inserted; `rng_type` and `value_type` must satisfy
     * seqan3::compatible_concept.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param count Number of copies.
     * \param value Element value to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `count==0`.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \par Complexity
     *
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     *
     * \par Example
     *
     * ```cpp
     * concatenated_sequences<dna4_vector> foobar;
     * foobar.insert(foobar.end(), 2, "ACGT"_dna4);
     * std::cout << foobar[0] << '\n'; // [A, C, G, T]
     * std::cout << foobar[1] << '\n'; // [A, C, G, T]
     * ```
     */
     template <typename rng_type>
    iterator insert(const_iterator pos, size_type const count, rng_type && value)
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)

    {
        auto const pos_as_num = std::distance(cbegin(), pos); // we want to insert BEFORE this position
        // TODO SEQAN_UNLIKELY
        if (count == 0)
            return begin() + pos_as_num;

        auto concatenated = ranges::view::repeat_n(value, count)
                          | ranges::view::join
                          | ranges::view::bounded;

        // we could call insert(pos, concatenated.cbegin(), concatenated.cend()) here
        // if view::join didn't strip the view of its random access and sized properties

        size_type value_len = 0;
        if constexpr (sized_range_concept<std::decay_t<rng_type>>)
            value_len = ranges::size(value);
        else
            value_len = std::distance(ranges::begin(value), ranges::end(value));

        data_values.reserve(data_values.size() + count * value_len);
        data_values.insert(data_values.cbegin() + data_delimiters[pos_as_num],
                           ranges::begin(concatenated),
                           ranges::end(concatenated));

        data_delimiters.reserve(data_values.size() + count);
        data_delimiters.insert(data_delimiters.cbegin() + pos_as_num,
                               count,
                               *(data_delimiters.cbegin() + pos_as_num));

        // adapt delimiters of inserted
        for (size_type i = 0; i < count; ++i)
            data_delimiters[pos_as_num + i + 1] += value_len * (i + 1);

        // adapt delimiters after that
        // TODO parallel execution policy or vectorization?
        std::for_each(data_delimiters.begin() + pos_as_num + count + 1,
                      data_delimiters.end(),
                      [full_len = value_len * count] (auto & d) { d += full_len; });

        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from range `[first, last)` before position in the container.
     * \tparam begin_iterator_type Must satisfy seqan3::forward_iterator_concept and satisfy seqan3::compatible_concept
     * with `value_type`.
     * \tparam end_iterator_type Must satisfy seqan3::sized_sentinel_concept.
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
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    template <typename begin_iterator_type, typename end_iterator_type>
    iterator insert(const_iterator pos, begin_iterator_type first, end_iterator_type last)
    //!\cond
        requires forward_iterator_concept<begin_iterator_type> &&
                 compatible_concept<begin_iterator_type, concatenated_sequences>
                 //&& sized_sentinel_concept<end_iterator_type, begin_iterator_type>
    //!\endcond
    {
        auto const pos_as_num = std::distance(cbegin(), pos);
        // TODO SEQAN_UNLIKELY
        if (last - first == 0)
            return begin() + pos_as_num;

        auto const ilist = ranges::make_iterator_range(first, last, std::distance(first, last));

        data_delimiters.reserve(data_values.size() + ilist.size());
        data_delimiters.insert(data_delimiters.cbegin() + pos_as_num,
                               ilist.size(),
                               *(data_delimiters.cbegin() + pos_as_num));


        // adapt delimiters of inserted region
        size_type full_len = 0;
        for (size_type i = 0; i < ilist.size(); ++i, ++first)
        {
            // constant for sized ranges and/or random access ranges, linear otherwise
            if constexpr (sized_range_concept<std::decay_t<decltype(*first)>>)
                full_len += ranges::size(*first);
            else
                full_len += std::distance(ranges::begin(*first), ranges::end(*first));

            data_delimiters[pos_as_num + 1 + i] += full_len;
        }

        // adapt values of inserted region
        auto concatenated = ilist | ranges::view::join | ranges::view::bounded;
        data_values.reserve(data_values.size() + full_len);
        data_values.insert(data_values.cbegin() + data_delimiters[pos_as_num],
                           ranges::begin(concatenated),
                           ranges::end(concatenated));

        // adapt delimiters behind inserted region
        // TODO parallel execution policy or vectorization?
        std::for_each(data_delimiters.begin() + pos_as_num + ilist.size() + 1,
                      data_delimiters.end(),
                      [full_len] (auto & d) { d += full_len; });

        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from initializer list before position in the container.
     * \tparam rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with `value_type`.
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
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    template <typename rng_type>
    iterator insert(const_iterator pos, std::initializer_list<rng_type> const & ilist)
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    {
        return insert(pos, ilist.begin(), ilist.end());
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
        auto const dist = std::distance(cbegin(), last);
        // TODO SEQAN_UNLIKELY
        if (last - first == 0)
            return begin() + dist;

        auto const distf = std::distance(cbegin(), first);

        // we need to scan once over the input
        size_type sum_size{0};
        for (; first != last; ++first)
            sum_size += ranges::size(*first);

        data_values.erase(data_values.cbegin() + data_delimiters[distf],
                          data_values.cbegin() + data_delimiters[dist]);

        data_delimiters.erase(data_delimiters.cbegin() + distf + 1,
                              data_delimiters.cbegin() + dist + 1);

        // adapt delimiters after that
        // TODO parallel execution policy or vectorization?
        std::for_each(data_delimiters.begin() + distf + 1,
                      data_delimiters.end(),
                      [sum_size] (auto & d) { d -= sum_size; });
        return begin() + dist;
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
       return erase(pos, pos + 1);
    }

    /*!\brief Appends the given element value to the end of the container.
     * \tparam rng_type The type of range to be inserted; must satisfy seqan3::compatible_concept with `value_type`.
     * \param value The value to append.
     *
     * If the new size() is greater than capacity() then all iterators and references (including the past-the-end
     * iterator) are invalidated. Otherwise only the past-the-end iterator is invalidated.
     *
     * \par Complexity
     *
     * Amortised linear in the size of value. Wort-case linear in concat_size().
     *
     * \par Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    template <typename rng_type>
    void push_back(rng_type && value)
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    {
        data_values.insert(data_values.cend(), ranges::begin(value), ranges::end(value));
        data_delimiters.push_back(data_delimiters.back() + ranges::size(value));
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
        assert(size() > 0);
        auto back_length = data_delimiters[size()] - data_delimiters[size() - 1];
        data_values.resize(data_values.size() - back_length);
        data_delimiters.pop_back();
    }

    /*!\brief Resizes the container to contain count elements.
     * \param count The new size.
     * \throws std::length_error If count > max_size().
     * \throws std::exception Any exception thrown by `Allocator::allocate()` (typically `std::bad_alloc`).
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
        assert(count < max_size());
        data_delimiters.resize(count + 1, data_delimiters.back());
        data_values.resize(data_delimiters.back());
    }

    /*!\copybrief resize()
     * \tparam rng_type The type of range to be inserted; `rng_type` and `value_type` must satisfy
     * seqan3::compatible_concept.
     * \param value Instead of appending empty containers, append copies of value.
     * \copydetails resize()
     */
    template <typename rng_type>
    void resize(size_type const count, rng_type && value)
        requires std::is_same_v<std::decay_t<rng_type>, std::decay_t<reference>> ||
                 std::is_same_v<std::decay_t<rng_type>, std::decay_t<const_reference>> ||
                 (forward_range_concept<std::decay_t<rng_type>> && compatible_concept<rng_type, value_type>)
    {
        assert(count < max_size());
        assert(concat_size() + count * ranges::size(value) < data_values.max_size());

        if (count < size())
            resize(count);
        else if (count > size())
            insert(cend(), count - size(), std::forward<rng_type>(value));
    }

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
    constexpr void swap(concatenated_sequences & rhs) noexcept
    {
        std::swap(data_values, rhs.data_values);
        std::swap(data_delimiters, rhs.data_delimiters);
    }

    //!\copydoc swap()
    constexpr void swap(concatenated_sequences && rhs) noexcept
    {
        std::swap(data_values, rhs.data_values);
        std::swap(data_delimiters, rhs.data_delimiters);
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(concatenated_sequences const & rhs) const noexcept
    {
        return data() == rhs.data();
    }

    constexpr bool operator!=(concatenated_sequences const & rhs) const noexcept
    {
        return data() != rhs.data();
    }

    constexpr bool operator<(concatenated_sequences const & rhs) const noexcept
    {
        return data() < rhs.data();
    }

    constexpr bool operator>(concatenated_sequences const & rhs) const noexcept
    {
        return data() > rhs.data();
    }

    constexpr bool operator<=(concatenated_sequences const & rhs) const noexcept
    {
        return data() <= rhs.data();
    }

    constexpr bool operator>=(concatenated_sequences const & rhs) const noexcept
    {
        return data() >= rhs.data();
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
        archive(data_values, data_delimiters);
    }
    //!\endcond
};

} // namespace seqan3
