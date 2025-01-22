// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::fast_istreambuf_iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <string>
#include <vector>

#include <seqan3/io/stream/detail/stream_buffer_exposer.hpp>

namespace seqan3::detail
{
/*!\brief Functionally the same as std::istreambuf_iterator, but faster.
 * \ingroup io_stream
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

    //!\brief Place to store a range of characters that overlaps stream buffer boundaries.
    std::string overflow_buffer{};

public:
    /*!\name Associated types
     * \{
     */
    using difference_type = ptrdiff_t;                 //!< Defaults to ptrdiff_t.
    using value_type = char_t;                         //!< The char type of the stream.
    using reference = char_t;                          //!< The char type of the stream.
    using pointer = void;                              //!< Has no pointer type.
    using iterator_category = std::input_iterator_tag; //!< Pure input iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fast_istreambuf_iterator() noexcept = default;                                             //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator const &) noexcept = default;             //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator &&) noexcept = default;                  //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator const &) noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator &&) noexcept = default;      //!< Defaulted.
    ~fast_istreambuf_iterator() noexcept = default;                                            //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit fast_istreambuf_iterator(std::basic_streambuf<char_t, traits_t> & ibuf) :
        stream_buf{reinterpret_cast<stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
    {
        assert(stream_buf != nullptr);

        if (stream_buf->gptr() == stream_buf->egptr()) // If current get area is empty,
            stream_buf->underflow();                   // ensure the stream buffer has content on construction.
    }
    //!\}

    //!\brief Cache until `raw_record.size() - 1` occurrences of `field_sep` followed by `record_end` were found.
    template <typename record_type>
        requires std::same_as<std::ranges::range_value_t<record_type>, std::string_view>
    void cache_record_into(char const record_end, char const field_sep, record_type & raw_record)
    {
        bool has_overflowed = false;
        size_t old_count = 0;
        char * data_begin = stream_buf->gptr(); // point into stream buffer by default
        size_t const number_of_fields = raw_record.size();
        size_t number_of_seen_fields = 0;
        std::vector<size_t> field_positions(number_of_fields, 0u);

        char const * ptr = stream_buf->gptr();

        auto overflow_into_buffer = [&]()
        {
            size_t count = stream_buf->egptr() - stream_buf->gptr();
            has_overflowed = true;
            overflow_buffer.resize(old_count + count);
            std::ranges::copy(stream_buf->gptr(), stream_buf->egptr(), overflow_buffer.data() + old_count);

            old_count += count;
            stream_buf->gbump(count);
            stream_buf->underflow();
        };

        while (number_of_seen_fields < number_of_fields - 1)
        {
            ptr = std::find(ptr, static_cast<char const *>(stream_buf->egptr()), field_sep);

            if (ptr != stream_buf->egptr()) // found an end of field
            {
                field_positions[number_of_seen_fields] = ptr - stream_buf->gptr() + old_count;
                ++ptr;
                ++number_of_seen_fields;
            }
            else
            {
                overflow_into_buffer();
                assert(stream_buf->gptr() != stream_buf->egptr()); // stream is not at end after overflow
                ptr = stream_buf->gptr();
            }
        }

        size_t count = 0;

        while (true) // Note: Might run idefinitely in release mode if no record_end is in input.
        {
            ptr = std::find(ptr, static_cast<char const *>(stream_buf->egptr()), record_end);

            if (ptr == stream_buf->egptr()) // stop_chr could not be found in current buffer
            {
                overflow_into_buffer();
                assert(stream_buf->gptr() != stream_buf->egptr()); // stream is not at end after overflow
                ptr = stream_buf->gptr();
            }
            else
            {
                count = ptr - stream_buf->gptr(); // processed characters until stop_chr has been found
                break;
            }
        }

        if (has_overflowed)
        {
            // need to copy last data
            overflow_buffer.resize(old_count + count);
            std::ranges::copy(stream_buf->gptr(), stream_buf->gptr() + count, overflow_buffer.data() + old_count);

            // make data pointer point into overflow
            data_begin = overflow_buffer.data();
        }

        stream_buf->gbump(count);

        // instantiate string_views in raw_record
        field_positions.back() = old_count + count;
        raw_record[0] = std::string_view{data_begin, field_positions[0]};
        for (size_t i = 1; i < number_of_fields; ++i)
            raw_record[i] = std::string_view{data_begin + field_positions[i - 1] + 1, data_begin + field_positions[i]};
    }

    //!\brief Cache `size` bytes from input stream.
    std::string_view cache_bytes(int32_t const size)
    {
        std::string_view result;

        if (stream_buf->egptr() - stream_buf->gptr() >= size)
        {
            result = std::string_view{stream_buf->gptr(), stream_buf->gptr() + size};
            stream_buf->gbump(size);
        }
        else
        {
            overflow_buffer.resize(size);

            int32_t remaining_bytes{size};

            while (stream_buf->egptr() - stream_buf->gptr() < remaining_bytes) // still not fully in buffer
            {
                std::ranges::copy(stream_buf->gptr(),
                                  stream_buf->egptr(),
                                  overflow_buffer.data() + size - remaining_bytes);
                size_t const number_of_copied_bytes = stream_buf->egptr() - stream_buf->gptr();
                remaining_bytes -= number_of_copied_bytes;
                stream_buf->gbump(number_of_copied_bytes);
                stream_buf->underflow();
                assert((remaining_bytes == 0 || stream_buf->egptr() != stream_buf->gptr())
                       && "I still need to read characters but my stream is at end.");
            }

            if (remaining_bytes != 0) // In stream_buf but not yet copied to overflow_buffer
            {
                std::ranges::copy(stream_buf->gptr(),
                                  stream_buf->gptr() + remaining_bytes,
                                  overflow_buffer.data() + size - remaining_bytes);

                stream_buf->gbump(remaining_bytes);
            }

            result = {overflow_buffer.begin(), overflow_buffer.end()};
        }

        return result;
    }

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
        assert(stream_buf->gptr() != stream_buf->egptr());
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
