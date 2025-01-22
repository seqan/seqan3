// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief A constexpr vector implementation with dynamic size at compile time.
 */

#pragma once

#include <array>
#include <type_traits>

#if SEQAN3_WITH_CEREAL
#    include <cereal/types/array.hpp>
#endif // SEQAN3_WITH_CEREAL

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/detail/integer_traits.hpp>
#include <seqan3/utility/views/repeat_n.hpp>

namespace seqan3
{

/*!\brief A constexpr vector implementation with dynamic size at compile time.
 * \implements seqan3::reservible_container
 * \implements seqan3::cerealisable
 * \ingroup utility_container
 * \tparam value_type_ The underlying value type stored in the vector.
 * \tparam capacity_   The capacity of the constexpr vector.
 *
 *
 * This implementation of a vector can be constructed, accessed and modified at compile time.
 * It has a fixed capacity but a dynamic size and provides all functionality of a sequence container. Note
 * that it also models a reservable sequence container but all associated member functions are no-op because the
 * capacity is fixed.
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
template <typename value_type_, size_t capacity_>
class small_vector
{
private:
    //!\brief Whether this class has noexpect members or not
    static constexpr bool is_noexcept = std::is_nothrow_copy_constructible_v<value_type_>;

public:
    /*!\name Associated types
     * \{
     */
    /*!\brief The value_type type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using value_type = value_type_;

    /*!\brief The reference type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using reference = value_type &;

    /*!\brief The const_reference type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using const_reference = value_type const &;

    /*!\brief The iterator type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using iterator = value_type *;

    /*!\brief The const_iterator type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using const_iterator = value_type const *;

    /*!\brief The difference_type type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using difference_type = ptrdiff_t;

    /*!\brief The size_type type.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    using size_type = detail::min_viable_uint_t<capacity_>;

    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr small_vector() noexcept = default;                                 //!< Defaulted.
    constexpr small_vector(small_vector const &) noexcept = default;             //!< Defaulted.
    constexpr small_vector(small_vector &&) noexcept = default;                  //!< Defaulted.
    constexpr small_vector & operator=(small_vector const &) noexcept = default; //!< Defaulted.
    constexpr small_vector & operator=(small_vector &&) noexcept = default;      //!< Defaulted.
    ~small_vector() noexcept = default;                                          //!< Defaulted.

    /*!\brief Construct from an (smaller or equally sized) array over the same value type.
     * \param[in] array The array to copy the data_ from.
     *
     * ### Complexity
     *
     * Linear in the size of `array`.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    explicit constexpr small_vector(std::array<value_type, capacity_> const & array) noexcept(is_noexcept) :
        data_{array},
        sz{capacity_}
    {}

    //!\cond
    template <size_t capacity2>
    explicit constexpr small_vector(std::array<value_type, capacity2> const & array) noexcept(is_noexcept) :
        sz{capacity2}
    {
        static_assert(capacity2 <= capacity_, "You can only initialize from array that has smaller or equal capacity.");
        std::ranges::copy(array, data_.begin());
    }
    //!\endcond

    /*!\brief Construct from a (smaller or equally sized) built in array over the same value type.
     * \param[in] array The array to copy the data_ from.
     *
     * ### Complexity
     *
     * Linear in `capacity_`.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t capacity2>
    explicit constexpr small_vector(value_type const (&array)[capacity2]) noexcept(is_noexcept) : sz{capacity2}
    {
        static_assert(capacity2 <= capacity_, "You can only initialize from array that has smaller or equal capacity.");
        std::ranges::copy(array, data_.begin());
    }

    /*!\brief Construct from a list of values of value_type.
     * \tparam other_value_type A parameter pack where each type is equal to value_type.
     * \param[in] args The values to construct from.
     * ### Complexity
     *
     * Linear in the number of elements.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename... other_value_type>
        requires (std::same_as<value_type, other_value_type> && ...)
    constexpr small_vector(other_value_type... args) noexcept(is_noexcept) :
        data_{args...},
        sz{sizeof...(other_value_type)}
    {
        static_assert(sizeof...(other_value_type) <= capacity_, "Value list must not exceed the capacity.");
    }

    /*!\brief Construct from two iterators.
     * \tparam begin_it_type Must model std::forward_iterator and `value_type` must be constructible from
     *                             the reference type of begin_it_type.
     * \tparam   end_it_type Must satisfy std::sentinel_for.
     * \param[in]         begin_it Begin of range to construct/assign from.
     * \param[in]           end_it End of range to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the distance between `begin_it` and `end_it`.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <std::forward_iterator begin_it_type, typename end_it_type>
        requires std::sentinel_for<end_it_type, begin_it_type>
              && std::constructible_from<value_type, std::iter_reference_t<begin_it_type>>
    constexpr small_vector(begin_it_type begin_it, end_it_type end_it) noexcept(is_noexcept) : small_vector{}
    {
        assign(begin_it, end_it);
    }

    /*!\brief Construct from a different range.
     * \tparam other_range_t The type of range to be inserted; must satisfy std::ranges::input_range and `value_type`
     *                       must be constructible from std::ranges::range_reference_t<other_range_t>.
     * \param[in]      range The sequences to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the size of `range`.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1. This is a non-standard C++ extension.}
     */
    template <std::ranges::input_range other_range_t>
        requires (!std::is_same_v<std::remove_cvref_t<other_range_t>, small_vector>)
              && std::constructible_from<value_type, std::ranges::range_reference_t<other_range_t>>
    explicit constexpr small_vector(other_range_t && range) noexcept(is_noexcept) :
        small_vector{std::ranges::begin(range), std::ranges::end(range)}
    {}

    /*!\brief Construct with `n` times `value`.
     * \param[in] n     Number of elements.
     * \param[in] value The initial value to be assigned.
     *
     * ### Complexity
     *
     * Linear in `n`.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr small_vector(size_type n, value_type value) noexcept(is_noexcept) : small_vector{}
    {
        assign(n, value);
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
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr small_vector & operator=(std::initializer_list<value_type> ilist) noexcept(is_noexcept)
    {
        assign(std::ranges::begin(ilist), std::ranges::end(ilist));
        return *this;
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
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void assign(std::initializer_list<value_type> ilist) noexcept(is_noexcept)
    {
        assign(std::ranges::begin(ilist), std::ranges::end(ilist));
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
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void assign(size_type const count, value_type const value) noexcept(is_noexcept)
    {
        clear();
        auto tmp = views::repeat_n(value, count);
        assign(std::ranges::begin(tmp), std::ranges::end(tmp));
    }

    /*!\brief Assign from a different range.
     * \tparam other_range_t The type of range to be inserted; must satisfy std::ranges::input_range and `value_type`
     *                       must be constructible from std::ranges::range_reference_t<other_range_t>.
     * \param[in]      range The sequences to construct/assign from.
     *
     * ### Complexity
     *
     * Linear in the size of `range`.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1. This is a non-standard C++ extension.}
     */
    template <std::ranges::input_range other_range_t>
        requires std::constructible_from<value_type, std::ranges::range_reference_t<other_range_t>>
    constexpr void assign(other_range_t && range) noexcept(is_noexcept)
    {
        assign(std::ranges::begin(range), std::ranges::end(range));
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
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <std::forward_iterator begin_it_type, typename end_it_type>
        requires std::sentinel_for<end_it_type, begin_it_type>
              && std::constructible_from<value_type, std::iter_reference_t<begin_it_type>>
    constexpr void assign(begin_it_type begin_it, end_it_type end_it) noexcept(is_noexcept)
    {
        clear();
        insert(cbegin(), begin_it, end_it);
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns the begin to the string.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr iterator begin() noexcept
    {
        return data_.data();
    }

    //!\copydoc seqan3::small_vector::begin()
    constexpr const_iterator begin() const noexcept
    {
        return data_.data();
    }

    //!\copydoc seqan3::small_vector::begin()
    constexpr const_iterator cbegin() const noexcept
    {
        return data_.data();
    }

    /*!\brief Returns iterator past the end of the vector.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr iterator end() noexcept
    {
        return data_.data() + sz;
    }

    //!\copydoc seqan3::small_vector::end()
    constexpr const_iterator end() const noexcept
    {
        return data_.data() + sz;
    }

    //!\copydoc seqan3::small_vector::end()
    constexpr const_iterator cend() const noexcept
    {
        return data_.data() + sz;
    }
    //!\}

    /*!\name Element access
     * \{
     */
    /*!\brief Return the i-th element.
     * \param[in] i Index of the element to retrieve.
     * \throws std::out_of_range If you access an element behind the last.
     * \returns A reference to the value at position `i`.
     *
     * \attention This is the only function of this class that is **not** constexpr because it might throw.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws std::out_of_range if `i >= size()`.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    reference at(size_type const i)
    {
        if (i >= size()) // [[unlikely]]
        {
            throw std::out_of_range{"Trying to access element behind the last in small_vector."};
        }
        return (*this)[i];
    }

    //!\copydoc seqan3::small_vector::at()
    const_reference at(size_type const i) const
    {
        if (i >= size()) // [[unlikely]]
        {
            throw std::out_of_range{"Trying to access element behind the last in small_vector."};
        }
        return (*this)[i];
    }

    /*!\brief Return the i-th element.
     * \param i The element to retrieve.
     * \returns A reference to the value at position `i`.
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
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr reference operator[](size_type const i) noexcept
    {
        assert(i < size());
        return data_[i];
    }

    //!\copydoc seqan3::small_vector::operator[]()
    constexpr const_reference operator[](size_type const i) const noexcept
    {
        assert(i < size());
        return data_[i];
    }

    /*!\brief Return the first element. Calling front on an empty container is undefined.
     * \returns A reference to the value at the first position.
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
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr reference front() noexcept
    {
        assert(size() > 0);
        return (*this)[0];
    }

    //!\copydoc seqan3::small_vector::front()
    constexpr const_reference front() const noexcept
    {
        assert(size() > 0);
        return (*this)[0];
    }

    /*!\brief Return the last element.
     * \returns A reference to the value at the last position.
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
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr reference back() noexcept
    {
        assert(size() > 0);
        return (*this)[size() - 1];
    }

    //!\copydoc seqan3::small_vector::back()
    constexpr const_reference back() const noexcept
    {
        assert(size() > 0);
        return (*this)[size() - 1];
    }

    /*!\brief Direct access to the underlying array.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr value_type * data() noexcept
    {
        return data_.data();
    }

    //!\copydoc seqan3::small_vector::data()
    constexpr value_type const * data() const noexcept
    {
        return data_.data();
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
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr bool empty() const noexcept
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
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr size_type size() const noexcept
    {
        return sz;
    }

    /*!\brief Returns the maximum number of elements the container is able to hold and resolves to `capacity_`.
     * \returns The number of elements in the container.
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
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr size_type max_size() const noexcept
    {
        return capacity_;
    }

    /*!\brief Returns the number of elements that the container is able to hold and resolves to `capacity_`.
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
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr size_type capacity() const noexcept
    {
        return capacity_;
    }

    /*!\brief Since the capacity is fixed on compile time, this is a no-op.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void reserve(size_type) const noexcept
    {
        // no-op
    }

    /*!\brief Since the capacity is fixed on compile time, this is a no-op.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void shrink_to_fit() const noexcept
    {
        // no-op
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
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void clear() noexcept
    {
        sz = 0;
    }

    /*!\brief Inserts value before position in the container.
     * \param   pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param value Element value to insert.
     * \returns     Iterator pointing to the inserted value.
     *
     * Inserting a value although the maximum capacity is reached is undefined behaviour.
     *
     * ### Complexity
     *
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr iterator insert(const_iterator pos, value_type const value) noexcept(is_noexcept)
    {
        return insert(pos, 1, value);
    }

    /*!\brief Inserts count copies of value before position in the container.
     * \param   pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param count Number of copies.
     * \param value Element value to insert.
     * \returns     Iterator pointing to the first element inserted, or `pos` if `count==0`.
     *
     * If `size()` + `count` > `capacity()` this function results in undefined behaviour.
     *
     * ### Complexity
     *
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr iterator insert(const_iterator pos, size_type const count, value_type const value) noexcept(is_noexcept)
    {
        auto tmp = views::repeat_n(value, count);
        return insert(pos, std::ranges::begin(tmp), std::ranges::end(tmp));
    }

    /*!\brief Inserts elements from range `[begin_it, end_it)` before position in the container.
     * \tparam begin_it_type Must satisfy std::forward_iterator and the `value_type` must be constructible from
     *                       the reference type of begin_it_type.
     * \tparam   end_it_type Must satisfy std::sentinel_for.
     * \param[in]        pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param[in]   begin_it Begin of range to construct/assign from.
     * \param[in]     end_it End of range to construct/assign from.
     * \returns              Iterator pointing to the first element inserted, or `pos` if `begin_it==end_it`.
     *
     * The behaviour is undefined if begin_it and end_it are iterators into `*this` or if, given the size `n` of the
     * range represented by [begin_t, end_it), `size()` + `n` > `capacity()`.
     *
     * ### Complexity
     *
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <std::forward_iterator begin_it_type, typename end_it_type>
        requires std::sentinel_for<end_it_type, begin_it_type>
              && std::constructible_from<value_type, std::iter_reference_t<begin_it_type>>
    constexpr iterator insert(const_iterator pos, begin_it_type begin_it, end_it_type end_it) noexcept(is_noexcept)
    {
        auto const pos_as_num = std::ranges::distance(cbegin(), pos);
        auto const length = std::ranges::distance(begin_it, end_it);

        assert(pos_as_num + length <= capacity());

        if (length == 0)
            return begin(); // nothing to insert

        for (size_type i = sz + length - 1; i > pos_as_num + length - 1; --i)
            data_[i] = data_[i - length];

        SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(-Wstringop-overflow)
        std::ranges::copy(begin_it, end_it, &data_[pos_as_num]);
        SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_STOP
        sz += length;
        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from initializer list before position in the container.
     * \param   pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param ilist Initializer list with values to insert.
     * \returns     Iterator pointing to the first element inserted, or `pos` if `ilist` is empty.
     *
     * Given the size `n` of `ilist`, this function results in undefined behaviour if `size()` + `n` > `capacity()`.
     *
     * ### Complexity
     *
     * Worst-case linear in size().
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr iterator insert(const_iterator pos, std::initializer_list<value_type> const & ilist) noexcept(is_noexcept)
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
     * The iterator begin_it does not need to be dereferenceable if begin_it==end_it: erasing an empty range is a no-op.
     *
     * ### Complexity
     *
     * Linear in size().
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr iterator erase(const_iterator begin_it, const_iterator end_it) noexcept
    {
        if (begin_it >= end_it) // [[unlikely]]
            return begin() + std::ranges::distance(cbegin(), end_it);

        size_type const length = std::ranges::distance(begin_it, end_it);
        auto out_it = begin() + std::ranges::distance(cbegin(), begin_it);

        while (end_it != cend())
            *(out_it++) = *(end_it++);

        sz -= length;
        return begin() + std::ranges::distance(cbegin(), begin_it);
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
     * No-throw guarantee.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr iterator erase(const_iterator pos) noexcept
    {
        return erase(pos, pos + 1);
    }

    /*!\brief Appends the given element value to the end of the container.
     * \param value The value to append.
     *
     * If the new size() is greater than capacity() this is undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void push_back(value_type const value) noexcept
    {
        assert(sz < capacity_);
        data_[sz] = value;
        ++sz;
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
     * No-throw guarantee.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void pop_back() noexcept
    {
        assert(sz > 0);
        --sz;
    }

    /*!\brief Resizes the container to contain count elements.
     * \param[in] count The new size.
     *
     * If count is greater than capacity this is undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void resize(size_type const count) noexcept
    {
        assert(count <= capacity_);
        sz = count;
    }

    /*!\copybrief seqan3::small_vector::resize()
     * \param value Append copies of value when resizing.
     * \copydetails seqan3::small_vector::resize()
     */
    constexpr void resize(size_type const count, value_type const value) noexcept
    {
        assert(count <= capacity_);
        for (size_t i = sz; i < count; ++i)
            data_[i] = value;
        sz = count;
    }

    /*!\brief Swap contents with another instance.
     * \param rhs The other instance to swap with.
     *
     * ### Complexity
     *
     * Linear in the size of both containers.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    constexpr void swap(small_vector & rhs) noexcept(is_noexcept)
    {
        auto tmp = *this;

        data_ = rhs.data_;
        sz = rhs.sz;

        rhs.data_ = tmp.data_;
        rhs.sz = tmp.sz;
    }

    //!\overload
    constexpr void swap(small_vector && rhs) noexcept(is_noexcept)
    {
        data_ = rhs.data_;
        sz = rhs.sz;
    }
    //!\}

    /*!\brief Swap contents with another instance.
     * \param lhs The first instance.
     * \param rhs The other instance to swap with.
     *
     * ### Complexity
     *
     * Linear in the size of both containers.
     *
     * ### Exceptions
     *
     * No-throw guarantee if value_type is std::is_nothrow_copy_constructible.
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    friend constexpr void swap(small_vector & lhs, small_vector & rhs) noexcept(is_noexcept)
    {
        lhs.swap(rhs);
    }

    //!\overload
    friend constexpr void swap(small_vector && lhs, small_vector && rhs) noexcept(is_noexcept)
    {
        lhs.swap(rhs);
    }

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Performs element-wise comparison.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t cap2>
    friend constexpr bool operator==(small_vector const & lhs, small_vector<value_type, cap2> const & rhs) noexcept
    {
        return std::ranges::equal(lhs, rhs);
    }

    /*!\brief Performs element-wise comparison.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t cap2>
    friend constexpr bool operator!=(small_vector const & lhs, small_vector<value_type, cap2> const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    /*!\brief Performs element-wise comparison.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t cap2>
    friend constexpr bool operator<(small_vector const & lhs, small_vector<value_type, cap2> const & rhs) noexcept
    {
        for (size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i)
            if (lhs[i] > rhs[i])
                return false;
            else if (lhs[i] < rhs[i])
                return true;
        return lhs.size() < rhs.size();
    }

    /*!\brief Performs element-wise comparison.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t cap2>
    friend constexpr bool operator>(small_vector const & lhs, small_vector<value_type, cap2> const & rhs) noexcept
    {
        for (size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i)
            if (lhs[i] < rhs[i])
                return false;
            else if (lhs[i] > rhs[i])
                return true;
        return lhs.size() > rhs.size();
    }

    /*!\brief Performs element-wise comparison.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t cap2>
    friend constexpr bool operator<=(small_vector const & lhs, small_vector<value_type, cap2> const & rhs) noexcept
    {
        return !(lhs > rhs);
    }

    /*!\brief Performs element-wise comparison.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <size_t cap2>
    friend constexpr bool operator>=(small_vector const & lhs, small_vector<value_type, cap2> const & rhs) noexcept
    {
        return !(lhs < rhs);
    }
    //!\}

public:
    //!\privatesection

    //!\brief Stores the actual data_.
    std::array<value_type, capacity_> data_{};
    //!\brief The size of the actual contained data_.
    size_type sz{0};

    //!\cond DEV
    /*!\brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(data_);
        archive(sz);
    }
    //!\endcond
};

/*!\name Type deduction guides
 * \{
 */

/*!\brief Deducts the size and value type from an built-in array on construction.
 * \details
 * \experimentalapi{Experimental since version 3.1.}
 */
template <size_t capacity2, typename value_type>
small_vector(value_type const (&array)[capacity2]) -> small_vector<value_type, capacity2>;
//!\}

} // namespace seqan3
