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

namespace seqan3::detail
{

/*!\brief The alignment
 * \ingroup execution
 * \tparam stream_buffer_t The buffer type for the alignment stream.
 *
 * \details
 *
 * Provides a stream-like range interface over the alignments instances that are computed in a
 * seqan3::detail::alignment_executor_two_way executor.
 */
template <typename stream_buffer_t>
class alignment_range
{
    static_assert(!std::is_const_v<stream_buffer_t>,
                  "Cannot create an alignment stream over a const buffer.");

    //!\brief The iterator of seqan3::detail::alignment_range.
    class iterator_type
    {
    public:
        //!\brief Type for distances between iterators.
        using difference_type = typename alignment_range::off_type;
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
        iterator_type()                                  = default;
        iterator_type(iterator_type const &)             = default;
        iterator_type(iterator_type &&)                  = default;
        iterator_type & operator=(iterator_type const &) = default;
        iterator_type & operator=(iterator_type &&)      = default;
        ~iterator_type()                                 = default;

        //!\brief Construct from alignment stream.
        iterator_type(alignment_range & _stream) : stream_ptr(&_stream)
        {}
        //!}

        /*!\name Read
         * \{
         */
        reference operator*()
        {
            return stream_ptr->cached();
        }

        const_reference operator*() const
        {
            return stream_ptr->cached();
        }
        //!\}

        /*!\name Increment operators
         * \{
         */
        iterator_type & operator++(/*pre*/)
        {
            stream_ptr->next();
            return *this;
        }

        iterator_type operator++(int /*post*/)
        {
            auto tmp{*this};
            stream_ptr->next();
            return tmp;
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */

        constexpr bool operator==(ranges::default_sentinel const & rhs) const
        {
            return stream_ptr->eof();
        }

        friend constexpr bool operator==(ranges::default_sentinel const & lhs,
                                         iterator_type const & rhs)
        {
            return rhs == lhs;
        }

        constexpr bool operator!=(ranges::default_sentinel const & rhs) const
        {
            return !(*this == rhs);
        }

        friend constexpr bool operator!=(ranges::default_sentinel const & lhs,
                                         iterator_type const & rhs)
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
    using off_type        = typename stream_buffer_t::off_type;
    //!\brief The alignment result type.
    using value_type      = typename stream_buffer_t::value_type;
    //!\brief The reference type.
    using reference       = typename stream_buffer_t::reference;
    //!\brief The iterator type.
    using iterator        = iterator_type;
    //!\brief The sentinel type.
    using sentinel        = ranges::default_sentinel;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_range()                                    = delete;
    alignment_range(alignment_range const &)             = default;
    alignment_range(alignment_range &&)                  = default;
    alignment_range & operator=(alignment_range const &) = default;
    alignment_range & operator=(alignment_range &&)      = default;
    ~alignment_range()                                   = default;

    alignment_range(stream_buffer_t & _stream_buffer) : executor_buffer{_stream_buffer}
    {}
    //!}

    /*!\name Iterators
     * \{
     */
    iterator begin()
    {
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
        if (auto opt = executor_buffer.bump(); opt.has_value())
            cache = &(*opt).get();
        else
            eof_flag = true;
    }

    //!\brief Returns the cached result.
    // TODO make iterable in case of multiple results per alignment.
    auto & cached() noexcept
    {
        assert(cache);
        return *cache;
    }

    //!\copydoc cached()
    auto const & cached() const noexcept
    {
        assert(cache);
        return *cache;
    }

    //!\brief Returns whether the executor buffer reached is end.
    constexpr bool eof() const noexcept
    {
        return eof_flag;
    }

private:
    //!\brief The underlying executor buffer.
    stream_buffer_t & executor_buffer;
    //!\brief Stores last read element.
    value_type      * cache{};
    //!\brief Indicates whether the stream has reached its end.
    bool eof_flag{false};
};

} // namespace seqan3
