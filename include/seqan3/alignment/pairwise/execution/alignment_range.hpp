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
 * \tparam alignment_buffer_t The buffer type for the alignment stream.
 *
 * \details
 *
 * Provides a stream-like range interface over the alignments instances that are computed in a
 * seqan3::detail::alignment_executor_two_way executor.
 */
template <typename alignment_buffer_t>
class alignment_range
{
    static_assert(!std::is_const_v<alignment_buffer_t>,
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
        //!\brief Use const reference type provided by container.
        using const_reference = std::add_const_t<reference>;
        //!\brief Pointer type is pointer of container element type.
        using pointer = std::add_pointer_t<value_type>;
        //!\brief Sets iterator category as input iterator.
        using iterator_category = std::input_iterator_tag;

        /*!\name Constructors, destructor and assignment
         * \{
         */
        iterator_type()                                  noexcept = default;
        iterator_type(iterator_type const &)             noexcept = default;
        iterator_type(iterator_type &&)                  noexcept = default;
        iterator_type & operator=(iterator_type const &) noexcept = default;
        iterator_type & operator=(iterator_type &&)      noexcept = default;
        ~iterator_type()                                          = default;

        //!\brief Construct from alignment stream.
        iterator_type(alignment_range & _stream) noexcept : stream_ptr(&_stream)
        {}
        //!}

        /*!\name Read
         * \{
         */
        reference operator*() noexcept
        {
            return stream_ptr->cache;
        }

        const_reference operator*() const noexcept
        {
            return stream_ptr->cache;
        }
        //!\}

        /*!\name Increment operators
         * \{
         */
        iterator_type & operator++(/*pre*/) noexcept
        {
            stream_ptr->next();
            return *this;
        }

        iterator_type operator++(int /*post*/) noexcept
        {
            auto tmp{*this};
            stream_ptr->next();
            return tmp;
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */

        constexpr bool operator==(std::ranges::default_sentinel const &) const noexcept
        {
            return stream_ptr->eof();
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
        alignment_range * stream_ptr{};
    };

    // Befriend the iterator with this class.
    friend class iterator_type;

public:

    //!\brief The offset type.
    using difference_type = typename alignment_buffer_t::difference_type;
    //!\brief The alignment result type.
    using value_type      = typename alignment_buffer_t::value_type;
    //!\brief The reference type.
    using reference       = typename alignment_buffer_t::reference;
    //!\brief The iterator type.
    using iterator        = iterator_type;
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

    //!\brief Explicit deletion to forbid copy construction of the underlying buffer.
    explicit alignment_range(alignment_buffer_t const & _alignment_buffer) = delete;

    /*!\brief Constructs a new alignment range by taking ownership over the passed alignment buffer.
     * \tparam _alignment_buffer_t   The buffer type. Must be the same type as `alignment_buffer_t`, when references
     *                               and cv-qualifiers are removed.
     * \param[in] _alignment_buffer  The buffer to take ownership from.
     *
     * \details
     *
     * Constructs a new alignment range by taking ownership over the passed alignment buffer.
     */
    explicit alignment_range(alignment_buffer_t && _alignment_buffer) :
        alignment_buffer{new alignment_buffer_t{std::move(_alignment_buffer)}},
        eof_flag(false)
    {}
    //!}

    /*!\name Iterators
     * \{
     */
    iterator begin()
    {
        if (!eof_flag)
            next();
        return iterator{*this};
    }

    sentinel end() noexcept
    {
        return {};
    }
    //!\}

protected:

    //!\brief Receives the next alignment result from the executor buffer.
    void next()
    {
        assert(!eof());

        if (!alignment_buffer)
            throw std::runtime_error{"No alignment execution buffer available."};

        if (auto opt = alignment_buffer->bump(); opt.has_value())
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
    std::unique_ptr<alignment_buffer_t> alignment_buffer{};
    //!\brief Stores last read element.
    value_type cache{};
    //!\brief Indicates whether the stream has reached its end.
    bool eof_flag{true};
};

/*!\name Type deduction guide
 * \relates seqan3::alignment_range
 * \{
 */

//!\brief Deduces from the passed alignment_buffer_t
template <typename alignment_buffer_t>
alignment_range(alignment_buffer_t &&) -> alignment_range<std::remove_reference_t<alignment_buffer_t>>;
//!}

} // namespace seqan3
