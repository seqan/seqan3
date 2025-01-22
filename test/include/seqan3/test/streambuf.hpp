// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides a stream class with a custom sized underlying buffer.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <streambuf>

#include <seqan3/core/platform.hpp>

namespace seqan3::test
{

/*!\brief Utility class to use a input streambuffer with custom buffer size.
 *
 * \attention Currently, this streambuf class only uses the small buffer for the get area (when using input streams).
 *            The class is extendable for output streams if needed.
 *
 * When workin with streams, it is sometimes necessary to test code with a small buffer s.t. edge cases are covered.
 * An edge case for SAM file IO is, for example, that a record spans the buffer boundary.
 *
 * ### input stream usage
 *
 * Use the class like this:
 *
 * ```
 * std::istringstream input{"This is what I want to read"};
 * std::istream & in{input};
 * std::streambuf * orig = in.rdbuf();
 * seqan3::test::small_buffer buf(orig);
 * in.rdbuf(&buf);
 *
 * // now use `in` whereever you would normally use `input` directly.
 * ```
 */
template <std::streamsize buffer_size>
class streambuf_with_custom_buffer_size : public std::streambuf
{
    std::streambuf * original_buffer;
    std::array<char, buffer_size> internal_buffer{};

public:
    streambuf_with_custom_buffer_size(std::streambuf * buf) : original_buffer(buf)
    {
        setg(internal_buffer.data(), internal_buffer.data() + buffer_size, internal_buffer.data() + buffer_size);
    }

    int underflow()
    {
        std::streamsize number_of_characters_read = original_buffer->sgetn(internal_buffer.data(), buffer_size);
        setg(internal_buffer.data(), internal_buffer.data(), internal_buffer.data() + number_of_characters_read);
        return (number_of_characters_read == buffer_size) ? *internal_buffer.data() : traits_type::eof();
    }
};

} // namespace seqan3::test
