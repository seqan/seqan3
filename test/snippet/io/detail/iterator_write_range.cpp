// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/stream/detail/fast_ostreambuf_iterator.hpp>

int main()
{
    std::string id{"seq1"};
    std::string sequence{"ACTGACTGACTGACTAGCATGACTAGCATGC"};

    // construct iterator from stream buffer
    std::ostringstream ostr{};
    auto stream_it = seqan3::detail::fast_ostreambuf_iterator{*ostr.rdbuf()};

    // You can do anything you could do with a regular std::ostreambuf_iterator
    stream_it = '>';  // writes '>' to stream
    *stream_it = ' '; // writes ' ' to stream

    // Additionally, there is an efficient write_range member function

    // Example 1: Write a range completely
    stream_it.write_range(id); // return value can be ignored

    // Example 2: Write a range in chunks of 10
    auto it = std::ranges::begin(sequence);
    while (it != std::ranges::end(sequence))
    {
        /* Note that you need cannot use stream_it.write_range(rng | std::views::take(10)) here
         * because the returned iterator is not of the correct type.
         */
        auto current_end = it;
        size_t steps = std::ranges::advance(current_end, 10u, std::ranges::end(sequence));
        using subrange_t =
            std::ranges::subrange<decltype(it), decltype(current_end), std::ranges::subrange_kind::sized>;
        // Be aware that your range_type must model std::ranges::borrowed_range in order to use the return value!
        it = stream_it.write_range(subrange_t{it, current_end, 10u - steps});
        stream_it = ' ';
    }
}
