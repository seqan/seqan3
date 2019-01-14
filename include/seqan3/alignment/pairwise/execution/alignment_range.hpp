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
 * \brief Provides seqan3::detail::alignment_range.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/utility/iterator.hpp>

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief The alignment
 * \ingroup execution
 * \tparam alignment_executor_type The buffer type for the alignment stream.
 *
 * \details
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
        constexpr iterator_type() noexcept = default;
        constexpr iterator_type(iterator_type const &) noexcept = default;
        constexpr iterator_type(iterator_type &&) noexcept = default;
        constexpr iterator_type & operator=(iterator_type const &) noexcept = default;
        constexpr iterator_type & operator=(iterator_type &&) noexcept = default;
        ~iterator_type() = default;

        //!\brief Construct from alignment stream.
        constexpr iterator_type(alignment_range & range) noexcept : range_ptr(&range)
        {}
        //!}

        /*!\name Read
         * \{
         */
        reference operator*() noexcept
        {
            return range_ptr->cache;
        }

        value_type const & operator*() const noexcept
        {
            return range_ptr->cache;
        }
        //!\}

        /*!\name Increment operators
         * \{
         */
        iterator_type & operator++(/*pre*/) noexcept
        {
            range_ptr->next();
            return *this;
        }

        void operator++(int /*post*/) noexcept
        {
            ++(*this);
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */

        constexpr bool operator==(std::ranges::default_sentinel const &) const noexcept
        {
            return range_ptr->eof();
        }

        friend constexpr bool operator==(std::ranges::default_sentinel const & lhs,
                                         iterator_type const & rhs) noexcept
        {
            return rhs == lhs;
        }

        constexpr bool operator!=(std::ranges::default_sentinel const & rhs) const noexcept
        {
            return !(*this == rhs);
        }

        friend constexpr bool operator!=(std::ranges::default_sentinel const & lhs,
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
    using sentinel        = std::ranges::default_sentinel;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_range()                                    = default;
    alignment_range(alignment_range const &)             = delete;
    alignment_range(alignment_range &&)                  = default;
    alignment_range & operator=(alignment_range const &) = delete;
    alignment_range & operator=(alignment_range &&)      = default;
    ~alignment_range()                                   = default;

    //!\brief Explicit deletion to forbid copy construction of the underlying executor.
    explicit alignment_range(alignment_executor_type const & _alignment_executor) = delete;

    /*!\brief Constructs a new alignment range by taking ownership over the passed alignment buffer.
     * \tparam _alignment_executor_type   The buffer type. Must be the same type as `alignment_executor_type`, when references
     *                               and cv-qualifiers are removed.
     * \param[in] _alignment_executor  The buffer to take ownership from.
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
    constexpr iterator begin()
    {
        if (!eof_flag)
            next();
        return iterator{*this};
    }

    const_iterator begin() const = delete;
    const_iterator cbegin() const = delete;

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
