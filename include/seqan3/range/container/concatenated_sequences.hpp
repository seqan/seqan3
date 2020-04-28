// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::concatenated_sequences.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>
#include <vector>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/views/as_const.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

#if SEQAN3_WITH_CEREAL
#include <cereal/types/vector.hpp>
#endif

namespace seqan3::detail
{

/*!\brief The reference type of seqan3::concatenated_sequences.
 * \tparam value_type The value_type of the seqan3::concatenated_sequences.
 * \tparam is_const_ref Reference type or const reference type.
 *
 * \details
 *
 * A light-weight type that inherits from the returned type of the seqan3::views::slice adaptor (std::span, std::ranges::subrange, ...),
 * but additionally provides implicit convertibility to the `value_type`. This is needed so that `value_type` and
 * `reference` type of seqan3::concatenated_sequences satisfy std::common_reference_with.
 *
 * The const version of this type additionally ensures deep constness to maintain container-like behaviour.
 */
template <typename value_type, bool const_>
struct concatenated_sequences_reference_proxy :
    public std::conditional_t<const_,
                              decltype(std::declval<value_type const &>() | views::as_const | views::slice(0,1)),
                              decltype(std::declval<value_type &>() | views::slice(0,1))>
{
    //!\brief The base type.
    using base_t =
        std::conditional_t<const_,
                           decltype(std::declval<value_type const &>() | views::as_const | views::slice(0,1)),
                           decltype(std::declval<value_type &>() | views::slice(0,1))>;

    //!\brief Inherit the base type's constructors.
    using base_t::base_t;

    //!\brief Construct from base type.
    concatenated_sequences_reference_proxy(base_t && rhs) : base_t{std::move(rhs)} {}

    //!\brief Implicitly convert to the `value_type` of seqan3::concatenated_sequences.
    operator value_type() const
    {
        value_type ret;
        ret.resize(std::ranges::size(*this));
        std::ranges::copy(*this, std::ranges::begin(ret));
        return ret;
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Container that stores sequences concatenated internally.
 * \tparam inner_type The type of sequences that will be stored. Must satisfy seqan3::reservible_container.
 * \tparam data_delimiters_type A container that stores the begin/end positions in the inner_type. Must be
 * seqan3::reservible_container and have inner_type's size_type as value_type.
 * \implements seqan3::cerealisable
 * \implements seqan3::reservible_container
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
 * * Also direct access to the delimiters via raw_data() [this is used by some algorithms].
 *
 * The disadvantages are:
 *
 * * Slower inserts and erases because the entire concatenation might have to be copied.
 * * No emplace operations.
 * * Modifying elements is limited to operations on elements of that element, i.e. you can change a character,
 * but you can't assign a new member sequence to an existing position.
 *
 * ###Example
 *
 * \include test/snippet/range/container/concatenated_sequences.cpp
 *
 * ###Exceptions
 *
 * Whenever a strong exception guarantee is given for this class, it presumes that
 * `std::is_nothrow_move_constructible<typename inner_type::value_type>` otherwise only basic exception safety can
 * be assumed.
 *
 * ###Thread safety
 *
 * This container provides no thread-safety beyond the promise given also by the STL that all
 * calls to `const` member function are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 */
template <typename inner_type,
          typename data_delimiters_type = std::vector<typename inner_type::size_type>>
//!\cond
    requires reservible_container<std::remove_reference_t<inner_type>> &&
             reservible_container<std::remove_reference_t<data_delimiters_type>> &&
             std::is_same_v<std::ranges::range_size_t<inner_type>, std::ranges::range_value_t<data_delimiters_type>>
//!\endcond
class concatenated_sequences
{
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

    //!\brief A proxy of type views::slice that represents the range on the concatenated vector.
    //!\hideinitializer
    using reference = detail::concatenated_sequences_reference_proxy<value_type, false>;

    //!\brief An immutable proxy of type views::slice that represents the range on the concatenated vector.
    //!\hideinitializer
    using const_reference = detail::concatenated_sequences_reference_proxy<value_type, true>;

    //!\brief The iterator type of this container (a random access iterator).
    //!\hideinitializer
    using iterator = detail::random_access_iterator<concatenated_sequences>;

    //!\brief The const iterator type of this container (a random access iterator).
    //!\hideinitializer
    using const_iterator = detail::random_access_iterator<concatenated_sequences const>;

    //!\brief A signed integer type (usually std::ptrdiff_t)
    //!\hideinitializer
    using difference_type = std::ranges::range_difference_t<data_delimiters_type>;

    //!\brief An unsigned integer type (usually std::size_t)
    //!\hideinitializer
    using size_type = std::ranges::range_size_t<data_delimiters_type>;
    //!\}

    //!\cond
    // this signals to range-v3 that something is a container :|
    using allocator_type    = void;
    //!\endcond

protected:
    /*!\name Compatibility
     * \brief Static constexpr variables that emulate/encapsulate seqan3::compatible (which doesn't work for types during their definition).
     * \{
     */
    //!\cond
    // unfortunately we cannot specialise the variable template so we have to add an auxiliary here
    template <std::ranges::range t>
    static constexpr bool is_compatible_with_value_type_aux(std::type_identity<t>)
    {
        return dimension_v<t> == dimension_v<value_type> &&
               std::convertible_to<std::ranges::range_reference_t<t>, std::ranges::range_value_t<value_type>>;
    }

    static constexpr bool is_compatible_with_value_type_aux(...)
    {
        return false;
    }
    //!\endcond

    //!\brief Whether a type satisfies seqan3::compatible with this class's `value_type` or `reference` type.
    //!\hideinitializer
    // we explicitly check same-ness, because these types may not be fully resolved, yet
    template <std::ranges::range t>
    static constexpr bool is_compatible_with_value_type = is_compatible_with_value_type_aux(std::type_identity<t>{});

    //!\brief Whether a type satisfies seqan3::compatible with this class.
    //!\hideinitializer
    // cannot use the concept, because this class is not yet fully defined
    template <typename t>
    //!\cond
        requires is_compatible_with_value_type<std::iter_reference_t<t>>
    //!\endcond
    static constexpr bool iter_value_t_is_compatible_with_value_type = true;

    //!\brief Whether a type satisfies seqan3::compatible with this class.
    //!\hideinitializer
    // cannot use the concept, because this class is not yet fully defined
    template <std::ranges::range t>
    //!\cond
        requires is_compatible_with_value_type<std::ranges::range_reference_t<t>>
    //!\endcond
    static constexpr bool range_value_t_is_compatible_with_value_type = true;
    //!\}

public:
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
     * \tparam rng_of_rng_type The type of range to be inserted; must satisfy
     *         \ref range_value_t_is_compatible_with_value_type.
     * \param rng_of_rng The sequences to construct/assign from.
     *
     * ###Complexity
     *
     * Linear in the cumulative size of `rng_of_rng`.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::input_range rng_of_rng_type>
    concatenated_sequences(rng_of_rng_type && rng_of_rng)
    //!\cond
        requires range_value_t_is_compatible_with_value_type<rng_of_rng_type>
    //!\endcond
    {
        if constexpr (std::ranges::sized_range<rng_of_rng_type>)
            data_delimiters.reserve(std::ranges::size(rng_of_rng) + 1);

        for (auto && val : rng_of_rng)
        {
            data_values.insert(data_values.end(), val.begin(), val.end());
            data_delimiters.push_back(data_delimiters.back() + val.size());
        }
    }

    /*!\brief Construct/assign with `count` times `value`.
     * \tparam rng_type The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param count Number of elements.
     * \param value The initial value to be assigned.
     *
     * ###Complexity
     *
     * In \f$O(count*value)\f$.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::forward_range rng_type>
    concatenated_sequences(size_type const count, rng_type && value)
    //!\cond
        requires is_compatible_with_value_type<rng_type>
    //!\endcond
    {
        // TODO SEQAN_UNLIKELY
        if (count == 0)
            return;

        insert(cend(), count, std::forward<rng_type>(value));
    }

    /*!\brief Construct/assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy std::forward_iterator and must satisfy
     *         \ref iter_value_t_is_compatible_with_value_type.
     * \tparam end_iterator_type Must satisfy std::sized_sentinel_for.
     * \param begin_it begin of range to construct/assign from.
     * \param end_it end of range to construct/assign from.
     *
     * ###Complexity
     *
     * Linear in the cumulative size of the ranges between `begin_it` and `end_it`.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::forward_iterator begin_iterator_type, typename end_iterator_type>
    concatenated_sequences(begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires std::sized_sentinel_for<end_iterator_type, begin_iterator_type> &&
                 iter_value_t_is_compatible_with_value_type<begin_iterator_type>
    //!\endcond
    {
        insert(cend(), begin_it, end_it);
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \tparam value_type_t The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param ilist an `std::initializer_list` of `value_type_t`.
     *
     * ###Complexity
     *
     * Linear in the cumulative size of the ranges in `ilist`.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::forward_range value_type_t = value_type>
    //!\cond
        requires is_compatible_with_value_type<value_type_t>
    //!\endcond
    concatenated_sequences(std::initializer_list<value_type_t> ilist)
    {
        assign(std::begin(ilist), std::end(ilist));
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \tparam value_type_t The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param ilist an `std::initializer_list` of `value_type_t`.
     *
     * ###Complexity
     *
     * Linear in the cumulative size of the elements in `ilist`.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::forward_range value_type_t>
    concatenated_sequences & operator=(std::initializer_list<value_type_t> ilist)
    //!\cond
        requires is_compatible_with_value_type<value_type_t>
    //!\endcond
    {
        assign(std::begin(ilist), std::end(ilist));
        return *this;
    }

    /*!\brief Construct/assign from a different range.
     * \tparam rng_of_rng_type The type of range to be inserted; must satisfy
     *         \ref range_value_t_is_compatible_with_value_type.
     * \param rng_of_rng The sequences to construct/assign from.
     *
     * ###Complexity
     *
     * Linear in the cumulative size of `rng_of_rng`.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::input_range rng_of_rng_type>
    void assign(rng_of_rng_type && rng_of_rng)
    //!\cond
        requires range_value_t_is_compatible_with_value_type<rng_of_rng_type>
    //!\endcond
    {
        concatenated_sequences rhs{std::forward<rng_of_rng_type>(rng_of_rng)};
        swap(rhs);
    }

    /*!\brief Construct/assign with `count` times `value`.
     * \tparam rng_type The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param count Number of elements.
     * \param value The initial value to be assigned.
     *
     * ###Complexity
     *
     * In \f$O(count*value)\f$.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::forward_range rng_type>
    void assign(size_type const count, rng_type && value)
    //!\cond
        requires (is_compatible_with_value_type<rng_type>)
    //!\endcond
    {
        concatenated_sequences rhs{count, value};
        swap(rhs);
    }

    /*!\brief Construct/assign from pair of iterators.
     * \tparam begin_iterator_type Must satisfy std::forward_iterator and satisfy
     *         \ref iter_value_t_is_compatible_with_value_type.
     * \tparam end_iterator_type Must satisfy std::sized_sentinel_for.
     * \param begin_it begin of range to construct/assign from.
     * \param end_it end of range to construct/assign from.
     *
     * ###Complexity
     *
     * Linear in the cumulative size of the ranges between `begin_it` and `end_it`.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::forward_iterator begin_iterator_type, typename end_iterator_type>
    void assign(begin_iterator_type begin_it, end_iterator_type end_it)
    //!\cond
        requires iter_value_t_is_compatible_with_value_type<begin_iterator_type> &&
                 std::sized_sentinel_for<end_iterator_type, begin_iterator_type>
    //!\endcond
    {
        concatenated_sequences rhs{begin_it, end_it};
        swap(rhs);
    }

    /*!\brief Construct/assign from `std::initializer_list`.
     * \tparam rng_type The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param ilist an `std::initializer_list` of `rng_type`.
     *
     * ###Complexity
     *
     * Linear in the cumulative size of the elements in `ilist`.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    template <std::ranges::forward_range rng_type = value_type>
    void assign(std::initializer_list<rng_type> ilist)
    //!\cond
        requires is_compatible_with_value_type<rng_type>
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (never modifies data)..
     */
    reference operator[](size_type const i)
    {
        assert(i < size());
        return data_values | views::slice(data_delimiters[i], data_delimiters[i+1]);
    }

    //!\copydoc operator[]()
    const_reference operator[](size_type const i) const
    {
        assert(i < size());
        return data_values | views::as_const | views::slice(data_delimiters[i], data_delimiters[i+1]);
    }

    /*!\brief Return the first element as a view. Calling front on an empty container is undefined.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Strong exception guarantee (never modifies data).
     */
    reference concat()
    {
        return data_values | views::slice(static_cast<size_type>(0), concat_size());
    }

    //!\copydoc concat()
    const_reference concat() const
    {
        return data_values | views::as_const | views::slice(static_cast<size_type>(0), concat_size());
    }

    /*!\brief Provides direct, unsafe access to underlying data structures.
     * \returns An std::pair of the concatenated sequences and the delimiter string.
     *
     * \details
     *
     * \noapi
     *
     * The exact representation of the data is implementation defined. Do not rely on it for API stability.
     */
    std::pair<decltype(data_values) &, decltype(data_delimiters) &> raw_data()
    {
        return {data_values, data_delimiters};
    }

    //!\copydoc raw_data()
    std::pair<decltype(data_values) const &, decltype(data_delimiters) const &> raw_data() const
    {
        return {std::as_const(data_values), std::as_const(data_delimiters)};
    }

    //!\copydoc raw_data()
    //!\deprecated Use raw_data() instead.
    SEQAN3_DEPRECATED_310 std::pair<decltype(data_values) &, decltype(data_delimiters) &> data()
    {
        return raw_data();
    }

    //!\copydoc raw_data()
    //!\deprecated Use raw_data() instead.
    SEQAN3_DEPRECATED_310 std::pair<decltype(data_values) const &, decltype(data_delimiters) const &> data() const
    {
        return raw_data();
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    /*!\brief Checks whether the container is empty.
     * \returns `true` if the container is empty, `false` otherwise.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * At most linear in the size() of the container.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * At most linear in the size() of the container.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * At most linear in the concat_size() of the container.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * \tparam rng_type The type of range to be inserted; must satisfy std::ranges::forward_range
     * and have the same `value_type` as `value_type` (i.e. `value_type`'s `value_type`!).
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param value Element value to insert.
     * \returns Iterator pointing to the inserted value.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ###Complexity
     *
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * ###Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     *
     * ###Example
     *
     * \include test/snippet/range/container/concatenated_sequences_insert.cpp
     */
    template <std::ranges::forward_range rng_type>
    iterator insert(const_iterator pos, rng_type && value)
    //!\cond
        requires is_compatible_with_value_type<rng_type>
    //!\endcond
    {
        return insert(pos, 1, std::forward<rng_type>(value));
    }
    // no specialisation for temporaries, since we have to copy anyway

    /*!\brief Inserts count copies of value before position in the container.
     * \tparam rng_type The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param count Number of copies.
     * \param value Element value to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `count==0`.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ###Complexity
     *
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * ###Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     *
     * ###Example
     *
     * \include test/snippet/range/container/concatenated_sequences_insert2.cpp
     */
    template <std::ranges::forward_range rng_type>
    iterator insert(const_iterator pos, size_type const count, rng_type && value)
    //!\cond
        requires is_compatible_with_value_type<rng_type>
    //!\endcond
    {
        auto const pos_as_num = std::distance(cbegin(), pos); // we want to insert BEFORE this position
        // TODO SEQAN_UNLIKELY
        if (count == 0)
            return begin() + pos_as_num;

        /* TODO implement views::flat_repeat_n that is like
         *  views::repeat_n(value, count) | views::join | ranges::view::bounded;
         * but preserves random access and size.
         *
         * then do
         *  auto concatenated = ranges::view::flat_repeat_n(value, count);
         *  insert(pos, concatenated.cbegin(), concatenated.cend())
         */

        size_type value_len = 0;
        if constexpr (std::ranges::sized_range<rng_type>)
            value_len = std::ranges::size(value);
        else
            value_len = std::distance(std::ranges::begin(value), std::ranges::end(value));

        data_values.reserve(data_values.size() + count * value_len);
        auto placeholder = views::repeat_n(std::ranges::range_value_t<rng_type>{}, count * value_len)
                         | std::views::common;
        // insert placeholder so the tail is moved once:
        data_values.insert(data_values.begin() + data_delimiters[pos_as_num],
                           std::ranges::begin(placeholder),
                           std::ranges::end(placeholder));

        // assign the actual values to the placeholder:
        size_t i = data_delimiters[pos_as_num];
        for (size_t j = 0; j < count; ++j)
            for (auto && v : value)
                data_values[i++] = v;

        data_delimiters.reserve(data_values.size() + count);
        data_delimiters.insert(data_delimiters.begin() + pos_as_num,
                               count,
                               *(data_delimiters.begin() + pos_as_num));

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
     * \tparam begin_iterator_type Must satisfy std::forward_iterator and
     *         \ref iter_value_t_is_compatible_with_value_type.
     * \tparam end_iterator_type Must satisfy std::sized_sentinel_for.
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
     * ###Complexity
     *
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * ###Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    template <std::forward_iterator begin_iterator_type, typename end_iterator_type>
    iterator insert(const_iterator pos, begin_iterator_type first, end_iterator_type last)
    //!\cond
        requires iter_value_t_is_compatible_with_value_type<begin_iterator_type> &&
                 std::sized_sentinel_for<end_iterator_type, begin_iterator_type>
    //!\endcond
    {
        auto const pos_as_num = std::distance(cbegin(), pos);
        // TODO SEQAN_UNLIKELY
        if (last - first == 0)
            return begin() + pos_as_num;

        auto const ilist =
            std::ranges::subrange<begin_iterator_type, end_iterator_type>(first,
                                                                          last,
                                                                          std::ranges::distance(first, last));

        data_delimiters.reserve(data_values.size() + ilist.size());
        data_delimiters.insert(data_delimiters.begin() + pos_as_num,
                               ilist.size(),
                               *(data_delimiters.begin() + pos_as_num));


        // adapt delimiters of inserted region
        size_type full_len = 0;
        for (size_type i = 0; i < ilist.size(); ++i, ++first)
        {
            full_len += std::ranges::distance(*first);
            data_delimiters[pos_as_num + 1 + i] += full_len;
        }

        // adapt values of inserted region
        auto placeholder = views::repeat_n(std::ranges::range_value_t<value_type>{}, full_len)
                         | std::views::common;
        // insert placeholder so the tail is moved only once:
        data_values.insert(data_values.begin() + data_delimiters[pos_as_num],
                           std::ranges::begin(placeholder),
                           std::ranges::end(placeholder));

        // assign the actual values to the placeholder:
        size_t i = data_delimiters[pos_as_num];
        for (auto && v0 : ilist)
            for (auto && v1 : v0)
                data_values[i++] = v1;


        // adapt delimiters behind inserted region
        // TODO parallel execution policy or vectorization?
        std::for_each(data_delimiters.begin() + pos_as_num + ilist.size() + 1,
                      data_delimiters.end(),
                      [full_len] (auto & d) { d += full_len; });

        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from initializer list before position in the container.
     * \tparam rng_type The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param ilist Initializer list with values to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `ilist` is empty.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * ###Complexity
     *
     * Worst-case linear in concat_size(). This is a drawback over e.g. `std::vector<std::vector<alphabet>>`.
     *
     * ###Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    template <std::ranges::forward_range rng_type>
    iterator insert(const_iterator pos, std::initializer_list<rng_type> const & ilist)
    //!\cond
        requires is_compatible_with_value_type<rng_type>
    //!\endcond
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
     * ###Complexity
     *
     * Linear in concat_size().
     *
     * ###Exceptions
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
            sum_size += std::ranges::size(*first);

        data_values.erase(data_values.begin() + data_delimiters[distf],
                          data_values.begin() + data_delimiters[dist]);

        data_delimiters.erase(data_delimiters.begin() + distf + 1,
                              data_delimiters.begin() + dist + 1);

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
     * ###Complexity
     *
     * Linear in concat_size().
     *
     * ###Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    iterator erase(const_iterator pos)
    {
       return erase(pos, pos + 1);
    }

    /*!\brief Appends the given element value to the end of the container.
     * \tparam rng_type The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param value The value to append.
     *
     * If the new size() is greater than capacity() then all iterators and references (including the past-the-end
     * iterator) are invalidated. Otherwise only the past-the-end iterator is invalidated.
     *
     * ###Complexity
     *
     * Amortised linear in the size of value. Wort-case linear in concat_size().
     *
     * ###Exceptions
     *
     * Basic exception guarantee, i.e. guaranteed not to leak, but container my contain invalid data after exceptions is
     * thrown.
     */
    template <std::ranges::forward_range rng_type>
    void push_back(rng_type && value)
    //!\cond
        requires is_compatible_with_value_type<rng_type>
    //!\endcond
    {
        data_values.insert(data_values.end(), std::ranges::begin(value), std::ranges::end(value));
        data_delimiters.push_back(data_delimiters.back() + std::ranges::size(value));
    }

    /*!\brief Removes the last element of the container.
     *
     * Calling pop_back on an empty container is undefined. In debug mode an assertion will be thrown.
     *
     * No iterators or references except for back() and end() are invalidated.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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
     * ###Complexity
     *
     * At most linear in the size() of the container.
     *
     * ###Exceptions
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
     * \tparam rng_type The type of range to be inserted; must satisfy \ref is_compatible_with_value_type.
     * \param value Instead of appending empty containers, append copies of value.
     * \copydetails resize()
     */
    template <std::ranges::forward_range rng_type>
    void resize(size_type const count, rng_type && value)
    //!\cond
        requires is_compatible_with_value_type<rng_type>
    //!\endcond
    {
        assert(count < max_size());
        assert(concat_size() + count * std::ranges::size(value) < data_values.max_size());

        if (count < size())
            resize(count);
        else if (count > size())
            insert(cend(), count - size(), std::forward<rng_type>(value));
    }

    /*!\brief Swap contents with another instance.
     * \param rhs The other instance to swap with.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
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

    /*!\name Comparison operators
     * \{
     */

    //!\brief Checks whether `*this` is equal to `rhs`.
    constexpr bool operator==(concatenated_sequences const & rhs) const noexcept
    {
        return raw_data() == rhs.raw_data();
    }

    //!\brief Checks whether `*this` is not equal to `rhs`.
    constexpr bool operator!=(concatenated_sequences const & rhs) const noexcept
    {
        return raw_data() != rhs.raw_data();
    }

    //!\brief Checks whether `*this` is less than `rhs`.
    constexpr bool operator<(concatenated_sequences const & rhs) const noexcept
    {
        return raw_data() < rhs.raw_data();
    }

    //!\brief Checks whether `*this` is greater than `rhs`.
    constexpr bool operator>(concatenated_sequences const & rhs) const noexcept
    {
        return raw_data() > rhs.raw_data();
    }

    //!\brief Checks whether `*this` is less than or equal to `rhs`.
    constexpr bool operator<=(concatenated_sequences const & rhs) const noexcept
    {
        return raw_data() <= rhs.raw_data();
    }

    //!\brief Checks whether `*this` is greater than or equal to `rhs`.
    constexpr bool operator>=(concatenated_sequences const & rhs) const noexcept
    {
        return raw_data() >= rhs.raw_data();
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
        archive(data_values, data_delimiters);
    }
    //!\endcond
};

} // namespace seqan3
