// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::fast_istreambuf_iterator and seqan3::fast_ostreambuf_iterator, as well as, seqan3::ostream
 *        and seqan3::ostreambuf iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iterator>

#ifndef __cpp_lib_ranges
#include <range/v3/iterator/stream_iterators.hpp>
#endif // __cpp_lib_ranges

#include <seqan3/core/platform.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
// ============================================================================
//  fast_istreambuf_iterator
// ============================================================================

/*!\brief Functionally the same as std::basic_streambuf<char_t, traits_t_>, but exposes protected members as public.
 * \ingroup stream
 * \tparam char_t   The stream's character type.
 * \tparam traits_t The stream's traits type.
 *
 * \details
 *
 * This wrapper adds no functionality to std::basic_streambuf and is only used to expose protected members to
 * access the get and put area of the std::basic_streambuf.
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
struct stream_buffer_exposer : public std::basic_streambuf<char_t, traits_t>
{
    //!\brief The actual stream type.
    using base_t = std::basic_streambuf<char_t, traits_t>;

    //!\cond
    // Expose protected members:
    using base_t::eback;
    using base_t::gptr;
    using base_t::egptr;
    using base_t::gbump;
    using base_t::underflow;

    using base_t::pbase;
    using base_t::pptr;
    using base_t::epptr;
    using base_t::pbump;
    using base_t::overflow;
    //!\endcond
};

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
    friend bool operator==(fast_istreambuf_iterator const & lhs, std::ranges::default_sentinel_t const &) noexcept
    {
        assert(lhs.stream_buf != nullptr);
        // compare size of remaining buffer; since ++ always resizes if possible, safe to compare pointers here
        return (lhs.stream_buf->gptr() == lhs.stream_buf->egptr());
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(fast_istreambuf_iterator const & lhs, std::ranges::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::ranges::default_sentinel);
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(std::ranges::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return rhs == std::ranges::default_sentinel;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(std::ranges::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return !(rhs == std::ranges::default_sentinel);
    }
    //!\}
};

/*!\brief Functionally the same as std::ostreambuf_iterator, but offers writing a range more efficiently.
 * \ingroup stream
 * \tparam char_t       The stream's character type.
 * \tparam traits_t     The stream's traits type.
 *
 * \details
 *
 * The functions seqan3::fast_ostreambuf_iterator::write_range and seqan3::fast_ostreambuf_iterator::write_n allow
 * more efficient writing of ranges by writing in chunks that avoiding overflow checks.
 *
 * \include test/snippet/io/detail/iterator_write_range.cpp
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
class fast_ostreambuf_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    stream_buffer_exposer<char_t, traits_t> * stream_buf = nullptr;

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;                //!< Defaults to ptrdiff_t.
    using value_type        = char_t;                   //!< The char type of the stream.
    using reference         = char_t;                   //!< The char type of the stream.
    using pointer           = void;                     //!< Has no pointer type.
    using iterator_category = std::output_iterator_tag; //!< Pure output iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fast_ostreambuf_iterator() noexcept = default; //!< Defaulted.
    fast_ostreambuf_iterator(fast_ostreambuf_iterator const &) noexcept = default; //!< Defaulted.
    fast_ostreambuf_iterator(fast_ostreambuf_iterator &&) noexcept = default; //!< Defaulted.
    fast_ostreambuf_iterator & operator=(fast_ostreambuf_iterator const &) noexcept = default; //!< Defaulted.
    fast_ostreambuf_iterator & operator=(fast_ostreambuf_iterator &&) noexcept = default; //!< Defaulted.
    ~fast_ostreambuf_iterator() noexcept = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit fast_ostreambuf_iterator(std::basic_streambuf<char_t, traits_t> & ibuf) :
        stream_buf{reinterpret_cast<stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
    {
        assert(stream_buf != nullptr);
        if (stream_buf->pptr() == stream_buf->epptr())
            stream_buf->overflow(); // ensures that put area has space available
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief no op.
    fast_ostreambuf_iterator & operator++()
    {
        return *this;
    }
    //!\overload
    fast_ostreambuf_iterator & operator++(int)
    {
        return *this;
    }
    //!\}

    //!\brief no op.
    fast_ostreambuf_iterator & operator*()
    {
        return *this;
    }

    //!\brief Writes a character to the associated output stream.
    fast_ostreambuf_iterator & operator=(char_t const c)
    {
        assert(stream_buf != nullptr);
        if (stream_buf->pptr() == stream_buf->epptr())
        {
            if (stream_buf->sputc(c) == traits_t::eof()) // overflow() [virtual], then write character
            {
                // LCOV_EXCL_START
                throw std::ios_base::failure{"Cannot write to output stream (reached traits::eof() condition)."};
                // LCOV_EXCL_STOP
            }
        }
        else
        {
            *stream_buf->pptr() = c;
            stream_buf->pbump(1); // advance pptr() in put area without any checks
        }
        return *this;
    }

    //!\brief Returns `true if this iterator has encountered the end-of-file condition on output, `false` otherwise.
    bool failed() const noexcept
    {
        return stream_buf->overflow() == traits_t::eof();
    }

    /*!\brief Writes a range to the associated output.
     * \tparam range_type The type of range to write; Must model std::ranges::forward_range.
     * \param[in] rng The range to write.
     * \returns If `range_type` models `std::ranges::forwarding_range` returns an iterator pointing to end of the range
     *          (rng) else returns `void`.
     *
     * This function avoids the buffer-at-end check by writing the range in chunks, where a chunks has the size of
     * the remaining space in the put area of the buffer.
     * If the range type models `std::ranges::sized_range` the chunks are written using `std::ranges::copy_n`, which
     * may use memcpy if applicable. Otherwise, a simple for loop iterates over the chunk.
     *
     * \attention You can only use the return value (end iterator) if your range type models
     *            `std::ranges::forwarding_range`.
     *
     * Example:
     *
     * \include test/snippet/io/detail/iterator_write_range.cpp
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires forwarding_range<range_type>
    //!\endcond
    auto write_range(range_type && rng)
    {
        using sen_t = std::ranges::sentinel_t<range_type>;
        using it_t = std::ranges::iterator_t<range_type>;

        it_t it = std::ranges::begin(rng);
        sen_t end = std::ranges::end(rng);

        while (it != end)
        {
            size_t const buffer_space = stream_buf->epptr() - stream_buf->pptr();

            if constexpr (std::ranges::sized_range<range_type>)
            {
                size_t const characters_to_write = std::min<size_t>(std::ranges::distance(it, end), buffer_space);
                auto copy_res = std::ranges::copy_n(it, characters_to_write, stream_buf->pptr());
                it = copy_res.in;
                stream_buf->pbump(characters_to_write);
            }
            else
            {
                size_t i = 0;
                for (; it != end && i < buffer_space; ++it, ++i)
                    *stream_buf->pptr() = *it;
                stream_buf->pbump(i);
            }

            if (it == end) // no more characters to write
                return it;

            // Push one more character and flush
            if (stream_buf->overflow(*it) == traits_t::eof())
            {
                // LCOV_EXCL_START
                throw std::ios_base::failure{"Cannot write to output stream (reached traits::eof() condition)."};
                // LCOV_EXCL_STOP
            }

            ++it; // drop 1 character that has been written in overflow()
        }

        return it;
    }

    //!\cond
    // overload for non-forwarding_range types that return void
    template <std::ranges::forward_range range_type>
    void write_range(range_type && rng)
    {
        write_range(rng); // lvalue is always a safe range. return value is ignored because iterator would be dangling
    }
    //!\endcond

    /*!\brief Writes a number to the underlying stream buffer using std::to_chars.
     * \tparam number_type The type of number; must model seqan3::arithmetic.
     * \param[in] num The number to write.
     */
    template <arithmetic number_type>
    auto write_number(number_type num)
    {
        if (stream_buf->epptr() - stream_buf->pptr() > 300) // enough space for any number, should be likely
        {
            auto res = std::to_chars(stream_buf->pptr(), stream_buf->epptr(), num);
            stream_buf->pbump(res.ptr - stream_buf->pptr()); // advance pptr
        }
        else
        {
            std::array<char, 300> arithmetic_buffer{};
            auto res = std::to_chars(&arithmetic_buffer[0], &arithmetic_buffer[0] + sizeof(arithmetic_buffer), num);
            write_range(std::ranges::subrange<char *, char *>(&arithmetic_buffer[0], res.ptr));
        }
    }

    /*!\brief Write `"\n"` or `"\r\n"` to the stream buffer, depending on arguments.
     * \param add_cr Whether to add carriage return, too.
     * \ingroup io
     */
    void write_end_of_line(bool const add_cr)
    {
        if (add_cr)
            *this = '\r';
        *this = '\n';
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\addtogroup stream
 * \{
 */

#ifdef __cpp_lib_ranges

/*!\typedef seqan3::ostream_iterator
* \brief Alias for std::ostream_iterator. Writes successive elements onto the output stream from which it was
*        constructed.
*/
using SEQAN3_DOXYGEN_ONLY(ostream_iterator =) ::std::ostream_iterator;

/*!\typedef seqan3::ostreambuf_iterator
 * \brief Alias for std::ostreambuf_iterator. Writes successive characters onto the output stream from which it
 *        was constructed.
 */
using SEQAN3_DOXYGEN_ONLY(ostreambuf_iterator =) ::std::ostreambuf_iterator;

#else

/*!\typedef seqan3::ostream_iterator
 * \brief Alias for ranges::ostream_iterator. Writes successive elements onto the output stream from which it was
 *        constructed.
 */
using SEQAN3_DOXYGEN_ONLY(ostream_iterator =) ::ranges::ostream_iterator;

/*!\typedef seqan3::ostreambuf_iterator
 * \brief Alias for ranges::ostreambuf_iterator. Writes successive characters onto the output stream from which it
 *        was constructed.
 */
using SEQAN3_DOXYGEN_ONLY(ostreambuf_iterator =) ::ranges::ostreambuf_iterator;

#endif // __cpp_lib_ranges
//!\}
} // namespace seqan3
