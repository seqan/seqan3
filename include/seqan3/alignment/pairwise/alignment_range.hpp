// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_range.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief The alignment
 * \ingroup pairwise_alignment
 * \tparam alignment_executor_type The buffer type for the alignment stream.
 *
 * \details
 * \attention This class considers moveable-only ranges.
 *
 * Provides a stream-like range interface over the alignments instances that are computed in a
 * seqan3::detail::alignment_executor_two_way executor.
 */
template <typename alignment_executor_type>
//TODO requires alignment_executor_concept<alignment_executor_type>
class alignment_range
{
    static_assert(!std::is_const_v<alignment_executor_type>,
                  "Cannot create an alignment stream over a const buffer.");

    //!\brief The iterator of seqan3::detail::alignment_range.
    class iterator_type
    {
    public:
        //!\brief Type for distances between iterators.
        using difference_type = typename alignment_range::difference_type;
        //!\brief Value type of container elements.
        using value_type = typename alignment_range::value_type;
        //!\brief Use reference type defined by container.
        using reference = typename alignment_range::reference;
        //!\brief Pointer type is pointer of container element type.
        using pointer = std::add_pointer_t<value_type>;
        //!\brief Sets iterator category as input iterator.
        using iterator_category = std::input_iterator_tag;

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr iterator_type() noexcept = default;                                  //!< Defaulted
        constexpr iterator_type(iterator_type const &) noexcept = default;             //!< Defaulted
        constexpr iterator_type(iterator_type &&) noexcept = default;                  //!< Defaulted
        constexpr iterator_type & operator=(iterator_type const &) noexcept = default; //!< Defaulted
        constexpr iterator_type & operator=(iterator_type &&) noexcept = default;      //!< Defaulted
        ~iterator_type() = default;                                                    //!< Defaulted

        //!\brief Construct from alignment stream.
        constexpr iterator_type(alignment_range & range) noexcept : range_ptr(&range)
        {}
        //!}

        /*!\name Read
         * \{
         */

         /*!\brief Access the pointed-to element.
          * \return A reference to the current element.
          */
        reference operator*() const noexcept
        {
            return range_ptr->cache;
        }
        //!\}

        /*!\name Increment operators
         * \{
         */

        //!\brief Increments the iterator by one.
        iterator_type & operator++(/*pre*/) noexcept
        {
            range_ptr->next();
            return *this;
        }

        //!\brief Returns an iterator incremented by one.
        void operator++(int /*post*/) noexcept
        {
            ++(*this);
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */

        //!\brief Checks whether `*this` is equal to the sentinel.
        constexpr bool operator==(std::ranges::default_sentinel_t const &) const noexcept
        {
            return range_ptr->eof();
        }

        //!\brief Checks whether `lhs` is equal to `rhs`.
        friend constexpr bool operator==(std::ranges::default_sentinel_t const & lhs,
                                         iterator_type const & rhs) noexcept
        {
            return rhs == lhs;
        }

        //!\brief Checks whether `*this` is not equal to the sentinel.
        constexpr bool operator!=(std::ranges::default_sentinel_t const & rhs) const noexcept
        {
            return !(*this == rhs);
        }

        //!\brief Checks whether `lhs` is not equal to `rhs`.
        friend constexpr bool operator!=(std::ranges::default_sentinel_t const & lhs,
                                         iterator_type const & rhs) noexcept
        {
            return rhs != lhs;
        }
        //!\}
    private:
        //!\brief Pointer to the underlying range.
        alignment_range * range_ptr{};
    };

    // Befriend the iterator with this class.
    // TODO Check if this is necessary.
    friend class iterator_type;

public:

    //!\brief The offset type.
    using difference_type = typename alignment_executor_type::difference_type;
    //!\brief The alignment result type.
    using value_type      = typename alignment_executor_type::value_type;
    //!\brief The reference type.
    using reference       = typename alignment_executor_type::reference;
    //!\brief The iterator type.
    using iterator        = iterator_type;
    //!\brief This range is never const-iterable. The const_iterator is always void.
    using const_iterator  = void;
    //!\brief The sentinel type.
    using sentinel        = std::ranges::default_sentinel_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_range() = default;                                   //!< Defaulted
    alignment_range(alignment_range const &) = delete;             //!< This is a move-only type.
    alignment_range(alignment_range &&) = default;                 //!< Defaulted
    alignment_range & operator=(alignment_range const &) = delete; //!< This is a move-only type.
    alignment_range & operator=(alignment_range &&) = default;     //!< Defaulted
    ~alignment_range() = default;                                  //!< Defaulted

    //!\brief Explicit deletion to forbid copy construction of the underlying executor.
    explicit alignment_range(alignment_executor_type const & _alignment_executor) = delete;

    /*!\brief Constructs a new alignment range by taking ownership over the passed alignment buffer.
     * \tparam _alignment_executor_type The buffer type. Must be the same type as `alignment_executor_type` when
     *                                  references and cv-qualifiers are removed.
     * \param[in] _alignment_executor   The buffer to take ownership from.
     *
     * \details
     *
     * Constructs a new alignment range by taking ownership over the passed alignment buffer.
     */
    explicit alignment_range(alignment_executor_type && _alignment_executor) :
        alignment_executor{new alignment_executor_type{std::move(_alignment_executor)}},
        eof_flag(false)
    {}
    //!}

    /*!\name Iterators
     * \{
     */

    /*!\brief Returns an iterator to the first element of the alignment range.
     * \return An iterator to the first element.
     *
     * \details
     * Invocation of this function will trigger the computation of the first alignment.
     */
    constexpr iterator begin()
    {
        if (!eof_flag)
            next();
        return iterator{*this};
    }

    const_iterator begin() const = delete;
    const_iterator cbegin() const = delete;

    /*!\brief Returns a sentinel signaling the end of the alignment range.
     * \return a sentinel.
     *
     * \details
     * The alignment range is an input range and the end is reached when the internal buffer over the alignment
     * results has signaled end-of-stream.
     */
    constexpr sentinel end() noexcept
    {
        return {};
    }

    constexpr sentinel end() const = delete;
    constexpr sentinel cend() const = delete;
    //!\}

protected:

    //!\brief Receives the next alignment result from the executor buffer.
    void next()
    {
        assert(!eof());

        if (!alignment_executor)
            throw std::runtime_error{"No alignment execution buffer available."};

        if (auto opt = alignment_executor->bump(); opt.has_value())
            cache = std::move(*opt);
        else
            eof_flag = true;
    }

    //!\brief Returns whether the executor buffer reached is end.
    constexpr bool eof() const noexcept
    {
        return eof_flag;
    }

private:
    //!\brief The underlying executor buffer.
    std::unique_ptr<alignment_executor_type> alignment_executor{};
    //!\brief Stores last read element.
    value_type cache{};
    //!\brief Indicates whether the stream has reached its end.
    bool eof_flag{true};
};

/*!\name Type deduction guide
 * \relates seqan3::alignment_range
 * \{
 */

//!\brief Deduces from the passed alignment_executor_type
template <typename alignment_executor_type>
alignment_range(alignment_executor_type &&) -> alignment_range<std::remove_reference_t<alignment_executor_type>>;
//!}

} // namespace seqan3
