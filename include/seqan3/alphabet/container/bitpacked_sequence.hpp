// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::bitpacked_sequence.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iterator>
#include <ranges>
#include <type_traits>

#include <sdsl/int_vector.hpp>

#include <seqan3/alphabet/detail/alphabet_proxy.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/to_rank.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/range/detail/random_access_iterator.hpp>
#include <seqan3/utility/math.hpp>
#include <seqan3/utility/views/convert.hpp>

namespace seqan3
{

/*!\brief A space-optimised version of std::vector that compresses multiple letters into a single byte.
 * \tparam alphabet_type The value type of the container, must satisfy seqan3::writable_semialphabet and std::regular.
 * \implements seqan3::reservible_container
 * \implements seqan3::cerealisable
 * \ingroup alphabet_container
 *
 * This class template behaves just like std::vector<alphabet_type> but has an internal representation where
 * multiple values are packed into a single byte/word to save space, e.g. seqan3::bitpacked_sequence<seqan3::dna4>
 * uses a quarter of of the memory that std::vector<seqan3::dna4> uses, because a single seqan3::dna4 letter can be
 * represented in two bits (instead of 8 which is the lower bound for a single object in C++).
 *
 * The disadvantages are slightly slower operations and unsafety towards parallel writes to adjacent positions
 * in the seqan3::bitpacked_sequence.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/container/bitpacked_sequence.cpp
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
 *
 * \stableapi{Since version 3.1.}
 */
template <writable_semialphabet alphabet_type>
    requires std::regular<alphabet_type>
class bitpacked_sequence
{
private:
    //!\brief The number of bits needed to represent a single letter of the alphabet_type.
    static constexpr size_t bits_per_letter = detail::ceil_log2(alphabet_size<alphabet_type>);

    static_assert(bits_per_letter <= 64, "alphabet must be representable in at most 64bit.");

    //!\brief Type of the underlying SDSL vector.
    using data_type = sdsl::int_vector<bits_per_letter>;

    //!\brief The data storage.
    data_type data;

    //!\brief Proxy data type returned by seqan3::bitpacked_sequence as reference to element.
    class reference_proxy_type : public alphabet_proxy<reference_proxy_type, alphabet_type>
    {
    private:
        //!\brief The base type.
        using base_t = alphabet_proxy<reference_proxy_type, alphabet_type>;
        //!\brief Befriend the base type so it can call our #on_update().
        friend base_t;

        //!\brief The proxy of the underlying data type.
        std::ranges::range_reference_t<data_type> internal_proxy;

        //!\brief Update the sdsl-proxy.
        constexpr void on_update() noexcept
        {
            internal_proxy = base_t::to_rank();
        }

    public:
        // Import from base:
        using base_t::operator=;

        /*!\name Constructors, destructor and assignment
         * \{
         */
        //!\brief Deleted, because using this proxy without a parent would be undefined behaviour.
        reference_proxy_type() = delete;
        constexpr reference_proxy_type(reference_proxy_type const &) noexcept = default;             //!< Defaulted.
        constexpr reference_proxy_type(reference_proxy_type &&) noexcept = default;                  //!< Defaulted.
        constexpr reference_proxy_type & operator=(reference_proxy_type const &) noexcept = default; //!< Defaulted.
        constexpr reference_proxy_type & operator=(reference_proxy_type &&) noexcept = default;      //!< Defaulted.
        ~reference_proxy_type() noexcept = default;                                                  //!< Defaulted.

        //!\brief Initialise from internal proxy type.
        reference_proxy_type(std::ranges::range_reference_t<data_type> const internal) noexcept :
            internal_proxy{internal}
        {
            // Call alphabet_base's assign_rank to prevent calling on_update() during construction
            // which is not necessary, because internal_proxy is already correctly initialised!
            base_t::base_t::assign_rank(static_cast<alphabet_rank_t<alphabet_type>>(internal));
        }
        //!\}
    };

    static_assert(writable_alphabet<reference_proxy_type>);
    //!\cond
    //NOTE(h-2): it is entirely unclear to me why we need this
    template <typename t>
        requires std::is_same_v<std::ranges::range_value_t<std::remove_cvref_t<t>>, alphabet_type>
    static constexpr bool has_same_value_type_v = true;
    //!\endcond

public:
    /*!\name Associated types
     * \{
     */
    /*!\brief Equals the alphabet_type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using value_type = alphabet_type;
    /*!\brief A proxy type (models seqan3::writable_semialphabet) that enables assignment, think of it as
     *        `value_type &`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using reference = reference_proxy_type;
    /*!\brief Equals the alphabet_type / value_type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using const_reference = alphabet_type;
    /*!\brief The iterator type of this container (a random access iterator).
     * \details
     * \stableapi{Since version 3.1.}
     */
    using iterator = detail::random_access_iterator<bitpacked_sequence>;
    /*!\brief The const_iterator type of this container (a random access iterator).
     * \details
     * \stableapi{Since version 3.1.}
     */
    using const_iterator = detail::random_access_iterator<bitpacked_sequence const>;
    /*!\brief A signed integer type (usually std::ptrdiff_t)
     * \details
     * \stableapi{Since version 3.1.}
     */
    using difference_type = std::ranges::range_difference_t<data_type>;
    /*!\brief An unsigned integer type (usually std::size_t)
     * \details
     * \stableapi{Since version 3.1.}
     */
    using size_type = std::ranges::range_size_t<data_type>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bitpacked_sequence() = default;                                                 //!< Defaulted.
    constexpr bitpacked_sequence(bitpacked_sequence const &) = default;             //!< Defaulted.
    constexpr bitpacked_sequence(bitpacked_sequence &&) = default;                  //!< Defaulted.
    constexpr bitpacked_sequence & operator=(bitpacked_sequence const &) = default; //!< Defaulted.
    constexpr bitpacked_sequence & operator=(bitpacked_sequence &&) = default;      //!< Defaulted.
    ~bitpacked_sequence() = default;                                                //!< Defaulted.

    /*!\brief Construct from a different range.
     * \tparam other_range_t The type of range to construct from; must satisfy std::ranges::input_range and
     *                       std::common_reference_with<std::ranges::range_value_t<other_range_t>, value_type>.
     * \param[in]      range The sequences to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the size of `range`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     *
     * \experimentalapi{Experimental since version 3.1. This is a non-standard C++ extension.}
     */
    template <typename other_range_t>
        requires (!std::same_as<bitpacked_sequence, std::remove_cvref_t<other_range_t>>)
              && std::ranges::input_range<other_range_t> && has_same_value_type_v<other_range_t>
    explicit bitpacked_sequence(other_range_t && range) :
        bitpacked_sequence{std::ranges::begin(range), std::ranges::end(range)}
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
     *
     * \stableapi{Since version 3.1.}
     */
    bitpacked_sequence(size_type const count, value_type const value) : data(count, to_rank(value))
    {}

    /*!\brief Construct from pair of iterators.
     * \tparam begin_iterator_type Must model std::forward_iterator and
     *                             std::common_reference_with<std::iter_value_t<begin_iterator_type>, value_type>.
     * \tparam   end_iterator_type Must model std::sentinel_for.
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
     *
     * \stableapi{Since version 3.1.}
     */
    template <std::forward_iterator begin_iterator_type, typename end_iterator_type>
    bitpacked_sequence(begin_iterator_type begin_it, end_iterator_type end_it)
        requires std::sentinel_for<end_iterator_type, begin_iterator_type>
              && std::common_reference_with<std::iter_value_t<begin_iterator_type>, value_type>
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
     *
     * \stableapi{Since version 3.1.}
     */
    bitpacked_sequence(std::initializer_list<value_type> ilist) : bitpacked_sequence(std::begin(ilist), std::end(ilist))
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
     *
     * \stableapi{Since version 3.1.}
     */
    bitpacked_sequence & operator=(std::initializer_list<value_type> ilist)
    {
        assign(std::begin(ilist), std::end(ilist));
        return *this;
    }

    /*!\brief Assign from a different range.
     * \tparam other_range_t The type of range to be inserted; must satisfy std::ranges::input_range and
     *                       std::common_reference_with<std::ranges::range_value_t<other_range_t>, value_type>.
     * \param[in]      range The sequences to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the size of `range`.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     *
     * \experimentalapi{Experimental since version 3.1. This is a non-standard C++ extension.}
     */
    template <std::ranges::input_range other_range_t>
    void assign(other_range_t && range)
        requires std::common_reference_with<std::ranges::range_value_t<other_range_t>, value_type>
    {
        bitpacked_sequence rhs{std::forward<other_range_t>(range)};
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
     *
     * \stableapi{Since version 3.1.}
     */
    void assign(size_type const count, value_type const value)
    {
        bitpacked_sequence rhs{count, value};
        swap(rhs);
    }

    /*!\brief Assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy std::forward_iterator and
     *                             std::common_reference_with<std::iter_value_t<begin_iterator_type>, value_type>.
     * \tparam   end_iterator_type Must satisfy std::sentinel_for.
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
     *
     * \stableapi{Since version 3.1.}
     */
    template <std::forward_iterator begin_iterator_type, typename end_iterator_type>
    void assign(begin_iterator_type begin_it, end_iterator_type end_it)
        requires std::sentinel_for<end_iterator_type, begin_iterator_type>
              && std::common_reference_with<std::iter_value_t<begin_iterator_type>, value_type>
    {
        bitpacked_sequence rhs{begin_it, end_it};
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
     */
    reference at(size_type const i)
    {
        if (i >= size()) // [[unlikely]]
        {
            throw std::out_of_range{"Trying to access element behind the last in bitpacked_sequence."};
        }
        return (*this)[i];
    }

    //!\copydoc at()
    const_reference at(size_type const i) const
    {
        if (i >= size()) // [[unlikely]]
        {
            throw std::out_of_range{"Trying to access element behind the last in bitpacked_sequence."};
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
     *
     * \stableapi{Since version 3.1.}
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
        return assign_rank_to(data[i], const_reference{});
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
     */
    reference back() noexcept
    {
        assert(size() > 0);
        return (*this)[size() - 1];
    }

    //!\copydoc back()
    const_reference back() const noexcept
    {
        assert(size() > 0);
        return (*this)[size() - 1];
    }

    /*!\brief Provides direct, unsafe access to underlying data structures.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * The exact representation of the data is implementation defined. Do not rely on it for API stability.
     *
     * \noapi
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}

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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
     */
    size_type max_size() const noexcept
    {
        return data.max_size();
    }

    /*!\brief Returns the number of elements that the container has currently allocated space for.
     * \returns The capacity of the currently allocated storage.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
     */
    void reserve(size_type const new_cap)
    {
        data.reserve(new_cap);
    }

    /*!\brief Requests the removal of unused capacity.
     *
     * It is a non-binding request to reduce capacity() to size().
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
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
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     *
     * \stableapi{Since version 3.1.}
     */
    iterator insert(const_iterator pos, size_type const count, value_type const value)
    {
        auto const pos_as_num = std::distance(cbegin(), pos); // we want to insert BEFORE this position

        data.insert(data.begin() + pos_as_num, count, to_rank(value));

        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from range `[begin_it, end_it)` before position in the container.
     * \tparam begin_iterator_type Must satisfy std::forward_iterator and
     *                             std::common_reference_with<std::iter_value_t<begin_iterator_type>, value_type>.
     * \tparam   end_iterator_type Must satisfy std::sentinel_for.
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
     *
     * \stableapi{Since version 3.1.}
     */
    template <std::forward_iterator begin_iterator_type, typename end_iterator_type>
    iterator insert(const_iterator pos, begin_iterator_type begin_it, end_iterator_type end_it)
        requires std::sentinel_for<end_iterator_type, begin_iterator_type>
              && std::common_reference_with<std::iter_value_t<begin_iterator_type>, value_type>
    {
        auto const pos_as_num = std::distance(cbegin(), pos);

        auto v = std::ranges::subrange<begin_iterator_type, end_iterator_type>{begin_it, end_it}
               | views::convert<value_type> | views::to_rank;
        data.insert(data.begin() + pos_as_num, std::ranges::begin(v), std::ranges::end(v));

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
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container may contain invalid data after exception is
     * thrown.
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
     */
    iterator erase(const_iterator begin_it, const_iterator end_it)
    {
        if (begin_it >= end_it) // [[unlikely]]
            return begin() + std::distance(cbegin(), end_it);

        auto const begin_it_pos = std::distance(cbegin(), begin_it);
        auto const end_it_pos = std::distance(cbegin(), end_it);

        data.erase(data.cbegin() + begin_it_pos, data.cbegin() + end_it_pos);

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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
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
     *
     * \stableapi{Since version 3.1.}
     */
    constexpr void swap(bitpacked_sequence & rhs) noexcept
    {
        std::swap(data, rhs.data);
    }

    //!\copydoc seqan3::bitpacked_sequence::swap
    constexpr void swap(bitpacked_sequence && rhs) noexcept
    {
        std::swap(data, rhs.data);
    }

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
     *
     * \stableapi{Since version 3.1.}
     */
    friend constexpr void swap(bitpacked_sequence & lhs, bitpacked_sequence & rhs) noexcept
    {
        std::swap(lhs, rhs);
    }

    //!\overload
    friend constexpr void swap(bitpacked_sequence && lhs, bitpacked_sequence && rhs) noexcept
    {
        std::swap(lhs, rhs);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    /*!\brief Checks whether `*this` is equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr bool operator==(bitpacked_sequence const & rhs) const noexcept
    {
        return data == rhs.data;
    }

    /*!\brief Checks whether `*this` is not equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr bool operator!=(bitpacked_sequence const & rhs) const noexcept
    {
        return data != rhs.data;
    }

    /*!\brief Checks whether `*this` is less than `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr bool operator<(bitpacked_sequence const & rhs) const noexcept
    {
        return data < rhs.data;
    }

    /*!\brief Checks whether `*this` is greater than `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr bool operator>(bitpacked_sequence const & rhs) const noexcept
    {
        return data > rhs.data;
    }

    /*!\brief Checks whether `*this` is less than or equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr bool operator<=(bitpacked_sequence const & rhs) const noexcept
    {
        return data <= rhs.data;
    }

    /*!\brief Checks whether `*this` is greater than or equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr bool operator>=(bitpacked_sequence const & rhs) const noexcept
    {
        return data >= rhs.data;
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(data);
    }
    //!\endcond
};

} // namespace seqan3
