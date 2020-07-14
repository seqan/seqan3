// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::fast_istreambuf_iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <seqan3/std/iterator>

#include <seqan3/io/stream/detail/stream_buffer_exposer.hpp>

namespace seqan3::detail
{
/*!\brief Functionally the same as std::istreambuf_iterator, but faster.
 * \ingroup stream
 * \tparam char_t       The stream's character type.
 * \tparam traits_t     The stream's traits type.
 *
 * \details
 *
 * Performs less virtual function calls than std::istreambuf_iterator.
 *
 * \todo Make this move-only after input iterators are allowed to be move-only.
 *
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
class fast_istreambuf_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    stream_buffer_exposer<char_t, traits_t> * stream_buf = nullptr;

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;               //!< Defaults to ptrdiff_t.
    using value_type        = char_t;                  //!< The char type of the stream.
    using reference         = char_t;                  //!< The char type of the stream.
    using pointer           = void;                    //!< Has no pointer type.
    using iterator_category = std::input_iterator_tag; //!< Pure input iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fast_istreambuf_iterator()                                              noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator const &)              noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator &&)                   noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator const &)  noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator &&)       noexcept = default; //!< Defaulted.
    ~fast_istreambuf_iterator()                                             noexcept = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit fast_istreambuf_iterator(std::basic_streambuf<char_t, traits_t> & ibuf) :
        stream_buf{reinterpret_cast<stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
    {
        assert(stream_buf != nullptr);
        stream_buf->underflow(); // ensure the stream buffer has content on construction
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advance by one and rebuffer if necessary (vtable lookup iff rebuffering).
    fast_istreambuf_iterator & operator++()
    {
        assert(stream_buf != nullptr);
        if ((stream_buf->gptr() + 1) == stream_buf->egptr())
            stream_buf->snextc(); // move right, then underflow()
        else
            stream_buf->gbump(1);
        return *this;
    }

    //!\overload
    void operator++(int)
    {
        ++(*this);
    }
    //!\}

    //!\brief Read current value from buffer (no vtable lookup, safe if not at end).
    reference operator*() const
    {
        assert(stream_buf != nullptr);
        return *stream_buf->gptr();
    }

    /*!\name Comparison operators
     * \brief We define comparison only against the sentinel.
     * \{
     */
    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(fast_istreambuf_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        assert(lhs.stream_buf != nullptr);
        // compare size of remaining buffer; since ++ always resizes if possible, safe to compare pointers here
        return (lhs.stream_buf->gptr() == lhs.stream_buf->egptr());
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(fast_istreambuf_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::default_sentinel);
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(std::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return rhs == std::default_sentinel;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(std::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return !(rhs == std::default_sentinel);
    }
    //!\}
};

} // namespace seqan3::detail
