// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::views::chunk.
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/range/concept.hpp>

namespace seqan3::detail
{
// ---------------------------------------------------------------------------------------------------------------------
// chunk_view class
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The type returned by seqan3::views::chunk.
 * \tparam urng_t The type of the underlying range, must model std::ranges::view and std::ranges::input_range.
 * \implements std::ranges::view
 * \implements std::ranges::random_access_range
 * \implements std::ranges::sized_range
 * \ingroup search_views
 *
 * \details
 *
 * Note that most members of this class are generated by std::ranges::view_interface which is not yet documented here.
 */
template <std::ranges::input_range urng_t>
    requires std::ranges::view<urng_t>
class chunk_view : public std::ranges::view_interface<chunk_view<urng_t>>
{
private:
    static_assert(std::ranges::input_range<urng_t>, "The chunk_view only works on input_ranges");

    //!\brief The underlying range.
    urng_t urange;

    //!\brief The chunk size to use.
    uint16_t chunk_size;

    // The iterator type if `urng_t` is a pure input range. See class definition for details.
    template <bool const_range>
    class basic_input_iterator;

    // The iterator type if `urng_t` is at least a forward range. See class definition for details.
    template <bool const_range>
    class basic_iterator;
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    chunk_view()
        requires std::default_initializable<urng_t>
    = default;                                                //!< Defaulted.
    chunk_view(chunk_view const & rhs) = default;             //!< Defaulted.
    chunk_view(chunk_view && rhs) = default;                  //!< Defaulted.
    chunk_view & operator=(chunk_view const & rhs) = default; //!< Defaulted.
    chunk_view & operator=(chunk_view && rhs) = default;      //!< Defaulted.
    ~chunk_view() = default;                                  //!< Defaulted.

    /*!\brief Construct from a view and the chunk size.
     * \param[in] underlying_range The underlying range to divide into chunks.
     * \param[in] size_of_chunk The size of the chunks, e.g. the length of the subrange returned at each position.
     */
    chunk_view(urng_t && underlying_range, uint16_t const size_of_chunk) : urange{std::forward<urng_t>(underlying_range)}, chunk_size{size_of_chunk}
    {}

    /*!\brief Construct from a non-view that can be view-wrapped and the chunk size.
     * \param[in] underlying_range The underlying range to divide into chunks.
     * \param[in] size_of_chunk The size of the chunks, e.g. the length of the subrange returned at each position.
     */
    template <typename rng_t>
    //!\cond
        requires (!std::same_as<std::remove_cvref_t<rng_t>, chunk_view>) && std::ranges::viewable_range<rng_t>
                  && std::constructible_from<urng_t, std::ranges::ref_view<std::remove_reference_t<rng_t>>>
    //!\endcond
    chunk_view(rng_t && underlying_range, uint16_t const size_of_chunk) :
        urange{std::views::all(std::forward<rng_t>(underlying_range))},
        chunk_size{size_of_chunk}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the range.
     * \returns Iterator to the first element.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant. For std::forward_ranges, O(c) when `c` is the chunk size, O(1) otherwise.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto begin() noexcept
    {
        if constexpr (std::ranges::forward_range<urng_t>)
            return basic_iterator<false>{std::ranges::begin(urange), std::ranges::end(urange), chunk_size};
        else
            return basic_input_iterator<false>{std::ranges::begin(urange), std::ranges::end(urange), chunk_size};
    }

    //!\copydoc begin()
    auto begin() const noexcept
        //!\cond
        requires const_iterable_range<urng_t>
    //!\endcond
    {
        if constexpr (std::ranges::forward_range<urng_t>)
            return basic_iterator<true>{std::ranges::cbegin(urange), std::ranges::cend(urange), chunk_size};
        else
            return basic_input_iterator<true>{std::ranges::cbegin(urange), std::ranges::cend(urange), chunk_size};
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * \details
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
    auto end() noexcept
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto end() const noexcept
        //!\cond
        requires const_iterable_range<urng_t>
    //!\endcond
    {
        return std::ranges::cend(urange);
    }
    //!\}

    /*!\brief Returns the size of the range, iff the underlying range is a std::ranges::sized_range.
     * \returns Size of range.
     */
    auto size()
        //!\cond
        requires std::ranges::sized_range<urng_t>
    //!\endcond
    {
        using size_type = std::ranges::range_size_t<urng_t>;
        return static_cast<size_type>((std::ranges::size(urange) + chunk_size - 1) / chunk_size); // round up
    }

    //!\copydoc size()
    auto size() const
        //!\cond
        requires std::ranges::sized_range<urng_t const>
    //!\endcond
    {
        using size_type = std::ranges::range_size_t<urng_t const>;
        return static_cast<size_type>((std::ranges::size(urange) + chunk_size - 1) / chunk_size); // round up
    }
};

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t>
chunk_view(rng_t &&, uint16_t const & chunk_size) -> chunk_view<std::views::all_t<rng_t>>;

// ---------------------------------------------------------------------------------------------------------------------
// chunk_view iterators (basic_input_iterator and basic_iterator)
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief Iterator for dividing an input range into chunks.
 * \tparam urng_t The type of the underlying range. Must model std::ranges::view and std::input_range but not std::forward_range.
 *
 * \details
 *
 * The basic_input_iterator can be used to iterate over an underlying input range in chunks.
 * It holds the start and end iterator of the underlying range, the chunk size and the subrange of the current chunk.
 * The current chunk is represented by a subrange on the underlying range with a std::views::take(chunk_size) applied.
 *
 * Depending on the type of the iterator passed to the basic_iterator, different functionality is available:
 *
 * | Concept modelled by passed text iterator | Available functions             |
 * |------------------------------------------|---------------------------------|
 * | std::input_iterator                      | \ref basic_input_iterator_comparison_chunk "Comparison operators"<br>\ref operator++ "Pre-increment (++it)"<br>\ref operator++(int) "Post-increment (it++)"<br>\ref operator* "Indirection operator (*it)" |
 *
 * \experimentalapi
 */
template <std::ranges::view urng_t>
template <bool const_range>
class chunk_view<urng_t>::basic_input_iterator :
    public maybe_iterator_category<maybe_const_iterator_t<const_range, urng_t>>
{
private:
    //!\brief The iterator type of the underlying range.
    using urng_it_t = maybe_const_iterator_t<const_range, urng_t>;

    //!\brief The sentinel type of the underlying range which is also the sentinel type of this iterator.
    using sentinel_t = maybe_const_sentinel_t<const_range, urng_t>;

    /*!\brief Helper iterator class to be used as iterator type in the subrange of this iterators value_type.
    * \tparam urng_t The type of the underlying range. Must model std::ranges::view and std::input_range but not std::forward_range.
    *
    * \details
    *
    * The only purpose of this class is to wrap the iterator type it is inheriting from (urng_it_t)
    * with the addition of decrementing the member variable `basic_input_iterator::remaining` when being incremented.
    * This way, the `basic_input_iterator` can keep track of how many times the `basic_input_iterator::urng_begin` has
    * been incremented within a chunk. A chunk can therefore be represented by a
    * `std::ranges::subrange<input_helper_iterator, sentinel_t>`.
    */
    template <typename outer_it_type>
    class input_helper_iterator : public urng_it_t
    {
    public:
        /*!\name Constructors, destructor and assignment
        * \{
        */
        constexpr input_helper_iterator() = default;                                          //!< Defaulted.
        constexpr input_helper_iterator(input_helper_iterator const &) = default;             //!< Defaulted.
        constexpr input_helper_iterator(input_helper_iterator &&) = default;                  //!< Defaulted.
        constexpr input_helper_iterator & operator=(input_helper_iterator const &) = default; //!< Defaulted.
        constexpr input_helper_iterator & operator=(input_helper_iterator &&) = default;      //!< Defaulted.
        ~input_helper_iterator() = default;                                                   //!< Defaulted.

        //!\brief Construct from the outer iterator and the underlying range iterator.
        input_helper_iterator(outer_it_type & outer_iterator, urng_it_t urng_it) : urng_it_t(urng_it)
        {
            outer_it = &outer_iterator;
        }

        //!\brief Construct from the underlying range iterator.
        input_helper_iterator(urng_it_t urng_it) : urng_it_t(urng_it)
        {}
        //!\}

        //!\brief Pre-increment will decrease the member variable basic_input_iterator::remaining.
        input_helper_iterator & operator++() noexcept
        {
            --(outer_it->remaining);
            urng_it_t::operator++();
            return *this;
        }

        //!\brief Post-increment will decrease the member variable basic_input_iterator::remaining.
        input_helper_iterator operator++(int) noexcept
        {
            input_helper_iterator tmp{*this};
            ++(*this);
            return tmp;
        }

        //!\brief Compare to the sentinel type (same as sentinel type of the underlying range).
        bool operator==(sentinel_t const & /* rhs */) noexcept
        {
            return this->outer_it->remaining == 0 || this->outer_it->urng_begin == this->outer_it->urng_end;
        }

        //!\brief Pointer to the outer iterator (basic_input_iterator).
        outer_it_type * outer_it{nullptr};
    };

    //!\brief This type will be used in the value_type of this iterator.
    using helper_it_t = input_helper_iterator<basic_input_iterator>;

    // befriend the iterator on a const range
    template <bool other_const_range>
    friend class basic_input_iterator;

    // befriend the inner iterator type
    template <typename outer_it_type>
    friend class input_helper_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = typename std::iter_difference_t<urng_it_t>;
    //!\brief Value type of this iterator.
    using value_type = std::ranges::subrange<helper_it_t, sentinel_t>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Same as `value_type`.
    using reference = value_type;
    //!\brief Tag this class is a pure input iterator.
    using iterator_concept = std::input_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr basic_input_iterator() = default;                                         //!< Defaulted.
    constexpr basic_input_iterator(basic_input_iterator const &) = default;             //!< Defaulted.
    constexpr basic_input_iterator(basic_input_iterator &&) = default;                  //!< Defaulted.
    constexpr basic_input_iterator & operator=(basic_input_iterator const &) = default; //!< Defaulted.
    constexpr basic_input_iterator & operator=(basic_input_iterator &&) = default;      //!< Defaulted.
    ~basic_input_iterator() = default;                                                  //!< Defaulted.

    //!\brief Allow iterator on a const range to be constructible from an iterator over a non-const range.
    constexpr basic_input_iterator(basic_input_iterator<!const_range> it) noexcept
        //!\cond
        requires const_range
    //!\endcond
    :
        chunk_size{std::move(it.chunk_size)},
        remaining{std::move(it.remaining)},
        urng_begin{std::move(it.urng_begin)},
        urng_end{std::move(it.urng_end)},
        current_chunk{std::move(it.current_chunk)}
    {}

    /*!\brief Construct from the start and end of the underlying range and a chunk size.
    * /param[in] it_begin Iterator pointing to the first position of the underlying range.
    * /param[in] it_end   Sentinel pointing to the end of the underlying range.
    * /param[in] size_of_chunk The chunk size, e.g. the length of the subrange returned by this iterator.
    *
    * \details
    *
    * ### Complexity
    *
    * Constant.
    */
    basic_input_iterator(urng_it_t it_begin, sentinel_t it_end, uint16_t size_of_chunk) :
        chunk_size{size_of_chunk},
        remaining{size_of_chunk},
        urng_begin{it_begin},
        urng_end{it_end}
    {
        current_chunk = std::ranges::subrange<helper_it_t, sentinel_t>{helper_it_t{*this, it_begin}, it_end};
    }
    //!\}

    //!\anchor basic_input_iterator_comparison_chunk
    //!\name Comparison operators
    //!\{

    //!\brief Compare to the sentinel type (same as sentinel type of the underlying range).
    friend bool operator==(basic_input_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return lhs.urng_begin == rhs;
    }

    //!\brief Compare to another basic_input_iterator.
    friend bool operator==(basic_input_iterator const & lhs, basic_input_iterator const & rhs) noexcept
    {
        return (lhs.urng_begin == rhs.urng_begin) && (lhs.remaining == rhs.remaining);
    }

    //!\brief Pre-increment.
    basic_input_iterator & operator++() noexcept
    {
        while (remaining > 0 && urng_begin != urng_end) // if chunk was not fully consumed and range is not at end
        {
            ++urng_begin;
            --remaining;
        }
        current_chunk = std::ranges::subrange<helper_it_t, sentinel_t>{helper_it_t{*this, urng_begin}, urng_end};
        remaining = chunk_size;
        return *this;
    }

    //!\brief Post-increment.
    basic_input_iterator operator++(int) noexcept
    {
        basic_input_iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    //!\brief Return the current chunk.
    value_type operator*() const noexcept
    {
        return current_chunk;
    }

private:
    //!\brief The chunk size, e.g. the length of the subrange returned by this iterator.
    uint16_t chunk_size;

    //!\brief The remaining elements in the chunk.
    uint16_t remaining;

    //!\brief Points to the start of the underlying range.
    urng_it_t urng_begin;

    //!\brief Points to the end of the underlying range.
    sentinel_t urng_end;

private:
    //!\brief The current chunk stored as a subrange.
    value_type current_chunk;
};

/*!\brief Iterator for dividing an input range into chunks.
 * \tparam urng_t The type of the underlying range. Must model std::ranges::view and std::forward_range.
 *
 * \details
 *
 * The basic_iterator can be used to iterate over an underlying input range in chunks.
 * It holds the start and end iterator of the underlying range, the chunk size and the subrange of the current chunk.
 * The current chunk is represented by a subrange on the underlying range.
 *
 * | Concept modelled by passed text iterator | Available functions             |
 * |------------------------------------------|---------------------------------|
 * | std::forward_iterator                    | \ref basic_iterator_comparison_chunk "Comparison operators"<br>\ref operator++ "Pre-increment (++it)"<br>\ref operator++(int) "Post-increment (it++)"<br>\ref operator* "Indirection operator (*it)" |
 * | std::bidirectional_iterator              | \ref operator-- "Pre-decrement (--it)"<br>\ref operator--(int) "Post-decrement (it--)" |
 * | std::random_access_iterator              | \ref operator+= "Forward (it +=)"<br>\ref operator+ "Forward copy (it +)"<br>\ref operator-= "Decrement(it -=)"<br>\ref basic_iterator_operator_decrement "Decrement copy (it -)"<br>\ref basic_iterator_operator-difference "Difference (it1 - it2)"<br>\ref operator[] "Subscript (it[])" |
 *
 * \experimentalapi
 */
template <std::ranges::view urng_t>
template <bool const_range>
class chunk_view<urng_t>::basic_iterator : public maybe_iterator_category<maybe_const_iterator_t<const_range, urng_t>>
{
private:
    //!\brief The iterator type of the underlying range.
    using it_t = maybe_const_iterator_t<const_range, urng_t>;
    //!\brief The sentinel type of the underlying range.
    using sentinel_t = maybe_const_sentinel_t<const_range, urng_t>;

    // befriend the iterator on a const range
    template <bool other_const_range>
    friend class basic_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = typename std::iter_difference_t<it_t>;
    //!\brief Value type of this iterator.
    using value_type = std::ranges::subrange<it_t, it_t>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Same as `value_type`.
    using reference = value_type;
    //!\brief Tag this class depending on which concept `it_t` models.
    using iterator_concept = std::conditional_t<std::contiguous_iterator<it_t>,
                                                typename std::random_access_iterator_tag,
                                                detail::iterator_concept_tag_t<it_t>>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr basic_iterator() = default;                                   //!< Defaulted.
    constexpr basic_iterator(basic_iterator const &) = default;             //!< Defaulted.
    constexpr basic_iterator(basic_iterator &&) = default;                  //!< Defaulted.
    constexpr basic_iterator & operator=(basic_iterator const &) = default; //!< Defaulted.
    constexpr basic_iterator & operator=(basic_iterator &&) = default;      //!< Defaulted.
    ~basic_iterator() = default;                                            //!< Defaulted.

    //!\brief Allow iterator on a const range to be constructible from an iterator over a non-const range.
    constexpr basic_iterator(basic_iterator<!const_range> const & it) noexcept
        //!\cond
        requires const_range
    //!\endcond
    :
        chunk_size{std::move(it.chunk_size)},
        urng_begin{std::move(it.urng_begin)},
        urng_end{std::move(it.urng_end)},
        current_chunk{std::move(it.current_chunk)}
    {}

    /*!\brief Construct from a given iterator on the text and a seqan3::shape.
    * /param[in] it_start Iterator pointing to the first position of the text.
    * /param[in] it_end   Sentinel pointing to the end of the text.
    * /param[in] size_of_chunk The chunk size, e.g. the length of the subrange returned by this iterator.
    *
    * \details
    *
    * ### Complexity
    *
    * Linear in chunk_size for pure forward ranges. Constant else.
    */
    basic_iterator(it_t it_start, sentinel_t it_end, uint16_t size_of_chunk) :
        chunk_size{size_of_chunk},
        urng_begin{it_start},
        urng_end{it_end}
    {
        current_chunk = value_type{it_start, get_next_end_of_chunk(it_start)};
    }
    //!\}

    //!\anchor basic_iterator_comparison_chunk
    //!\name Comparison operators
    //!\{

    //!\brief Compare to end of underlying range.
    friend bool operator==(basic_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return lhs.current_chunk.begin() == rhs;
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return (lhs.current_chunk.begin() == rhs.current_chunk.begin()) && (lhs.chunk_size == rhs.chunk_size);
    }

    //!\brief Compare to underlying range sentinal type.
    friend bool operator!=(basic_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator!=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator<(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return (lhs.chunk_size <= rhs.chunk_size) && (lhs.current_chunk.begin() < rhs.current_chunk.begin());
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator>(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return (lhs.chunk_size >= rhs.chunk_size) && (lhs.current_chunk.begin() > rhs.current_chunk.begin());
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator<=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return (lhs.chunk_size <= rhs.chunk_size) && (lhs.current_chunk.begin() <= rhs.current_chunk.begin());
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator>=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return (lhs.chunk_size >= rhs.chunk_size) && (lhs.current_chunk.begin() >= rhs.current_chunk.begin());
    }

    //!\}

    //!\brief Pre-increment.
    basic_iterator & operator++() noexcept
    {
        current_chunk = value_type{current_chunk.end(), get_next_end_of_chunk(current_chunk.end())};
        return *this;
    }

    //!\brief Post-increment.
    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    /*!\brief Pre-decrement.
     * \attention This function is only available if `it_t` models std::bidirectional_iterator.
     */
    basic_iterator & operator--() noexcept
        //!\cond
        requires std::bidirectional_iterator<it_t>
    //!\endcond
    {
        current_chunk = value_type{get_former_start_of_chunk(current_chunk.begin()), current_chunk.begin()};
        return *this;
    }

    /*!\brief Post-decrement.
     * \attention This function is only available if `it_t` models std::bidirectional_iterator.
     */
    basic_iterator operator--(int) noexcept
        //!\cond
        requires std::bidirectional_iterator<it_t>
    //!\endcond
    {
        basic_iterator tmp{*this};
        --(*this);
        return tmp;
    }

    /*!\brief Forward this iterator.
     * \attention This function is only available if `it_t` models std::random_access_iterator.
     */
    basic_iterator & operator+=(difference_type const skip) noexcept
        //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        auto new_start_it = current_chunk.begin() + (chunk_size * skip);
        current_chunk = value_type{new_start_it, get_next_end_of_chunk(new_start_it)};
        return *this;
    }

    /*!\brief Forward copy of this iterator.
     * \attention This function is only available if `it_t` models std::random_access_iterator.
     */
    basic_iterator operator+(difference_type const skip) const noexcept
        //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        basic_iterator tmp{*this};
        return tmp += skip;
    }

    /*!\brief Non-member operator+ delegates to non-friend operator+.
     * \attention This function is only available if `it_t` models std::random_access_iterator.
     */
    friend basic_iterator operator+(difference_type const skip, basic_iterator const & it) noexcept
        //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        return it + skip;
    }

    /*!\brief Decrement iterator by `skip`.
     * \attention This function is only available if `it_t` models std::random_access_iterator.
     */
    basic_iterator & operator-=(difference_type const skip) noexcept
        //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        auto new_start_it = current_chunk.begin() - (chunk_size * skip);
        current_chunk = value_type{new_start_it, get_next_end_of_chunk(new_start_it)};
        return *this;
    }

    /*!\anchor basic_iterator_operator_decrement
     * \brief Return decremented copy of this iterator.
     * \attention This function is only available if `it_t` models std::random_access_iterator.
     */
    basic_iterator operator-(difference_type const skip) const noexcept
        //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        basic_iterator tmp{*this};
        return tmp -= skip;
    }

    /*!\brief Non-member operator- delegates to non-friend operator-.
     * \attention This function is only available if `it_t` models std::random_access_iterator.
     */
    friend basic_iterator operator-(difference_type const skip, basic_iterator const & it) noexcept
        //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        return it - skip;
    }

    /*!\anchor basic_iterator operator-difference
     * \brief Return offset between two iterator's positions.
     * \attention This function is only available if `it_t` models std::sized_sentinel_for<it_t, it_t>.
     */
    friend difference_type operator-(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
        //!\cond
        requires std::sized_sentinel_for<it_t, it_t>
    //!\endcond
    {
        return static_cast<difference_type>((lhs.current_chunk.begin() - rhs.current_chunk.begin()) / lhs.chunk_size);
    }

    /*!\brief Return offset between remote sentinel's position and this.
     * \attention This function is only available if sentinel_t and it_t model std::sized_sentinel_for.
     */
    friend difference_type operator-(sentinel_t const & /* lhs */, basic_iterator const & rhs) noexcept
        //!\cond
        requires std::sized_sentinel_for<sentinel_t, it_t>
    //!\endcond
    {
        return static_cast<difference_type>((rhs.urng_end - rhs.current_chunk.begin() + rhs.chunk_size - 1)
                                            / rhs.chunk_size);
    }

    /*!\brief Return offset this and remote sentinel's position.
     * \attention This function is only available if it_t and sentinel_t model std::sized_sentinel_for.
     */
    friend difference_type operator-(basic_iterator const & lhs, sentinel_t const & rhs) noexcept
        //!\cond
        requires std::sized_sentinel_for<sentinel_t, it_t>
    //!\endcond
    {
        return -(rhs - lhs);
    }

    /*!\brief Move the iterator by a given offset and return the corresponding chunk (subrange).
     * \attention This function is only available if `it_t` models std::random_access_iterator.
     */
    reference operator[](difference_type const n) const
        //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        return *(*this + n);
    }

    //!\brief Return the current chunk, e.g the current subrange.
    value_type operator*() const noexcept
    {
        return current_chunk;
    }

private:
    //!\brief The chunk size, e.g. the length of the subrange returned by this iterator.
    uint16_t chunk_size;

    //!\brief Points to the start of the underlying range.
    it_t urng_begin;

    //!\brief Points to the end of the underlying range.
    sentinel_t urng_end;

    //!\brief The current chunk stored as a subrange.
    value_type current_chunk;

    //!\brief Move to the end of the next chunk.
    it_t get_next_end_of_chunk(it_t start_of_chunk)
    {
        it_t end_of_chunk{start_of_chunk};
        for (uint16_t i = 0; i < chunk_size && end_of_chunk != urng_end; ++i)
            ++end_of_chunk;
        return end_of_chunk;
    }

    //!\brief Move to the start of the former chunk.
    it_t get_former_start_of_chunk(it_t end_of_chunk)
    {
        it_t start_of_chunk{end_of_chunk};
        for (uint16_t i = 0; i < chunk_size && start_of_chunk != urng_begin; ++i)
            --start_of_chunk;
        return start_of_chunk;
    }
};

// ---------------------------------------------------------------------------------------------------------------------
// chunk_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//!\brief views::chunk's range adaptor object type (non-closure).
//!\ingroup search_views
struct chunk_fn
{
    //!\brief Store the `chunk_size` and return a range adaptor closure object.
    constexpr auto operator()(uint16_t const chunk_size) const
    {
        return adaptor_from_functor{*this, chunk_size};
    }

    /*!\brief Call the view's constructor with the underlying range and a chunk_size as argument.
     * \param[in] urange The range to process. Must model std::ranges::range.
     * \param[in] chunk_size TThe chunk size, e.g. the length of the subrange returned by this iterator.
     * \returns A range of subranges.
     */
    template <std::ranges::range underlying_range_t>
    constexpr auto operator()(underlying_range_t && urange, uint16_t const chunk_size) const
    {
        static_assert(std::ranges::viewable_range<underlying_range_t>,
                      "The range parameter to views::chunk cannot be a temporary of a non-view range.");
        static_assert(std::ranges::input_range<underlying_range_t>,
                      "The range parameter to views::chunk must model std::ranges::input_range.");

        return chunk_view{std::forward<underlying_range_t>(urange), chunk_size};
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{
/*!\brief                Divide a range in chunks.
 * \tparam urng_t        The type of the range being processed. See below for requirements. [template parameter is
 *                       omitted in pipe notation]
 * \param[in] urange     The range being processed. [parameter is omitted in pipe notation]
 * \param[in] chunk_size The seqan3::shape that determines how to compute the hash value.
 * \returns              A range of subranges pointing to the underlying range.
 *                       See below for the properties of the returned range.
 * \ingroup utility_views
 *
 * \details
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       |                                    | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *preserved*                      |
 * | std::ranges::random_access_range |                                    | *preserved*                      |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *preserved*                      |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   |                                    | todo                             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/utility/views/chunk.cpp
 *
 * \hideinitializer
 *
 * \experimentalapi
 */
inline constexpr auto chunk = detail::chunk_fn{};

} // namespace seqan3::views
