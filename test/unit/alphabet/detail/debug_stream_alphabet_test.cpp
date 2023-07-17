// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>

struct gap_alphabet
{
    using rank_t = uint8_t;
    using char_t = char;

    rank_t to_rank() const noexcept
    {
        return 0;
    }
    char_t to_char() const noexcept
    {
        return '-';
    }

    static constexpr bool alphabet_size{1};

    friend bool operator<(gap_alphabet, gap_alphabet)
    {
        return false;
    }
    friend bool operator<=(gap_alphabet, gap_alphabet)
    {
        return true;
    }
    friend bool operator>(gap_alphabet, gap_alphabet)
    {
        return false;
    }
    friend bool operator>=(gap_alphabet, gap_alphabet)
    {
        return true;
    }
    friend bool operator==(gap_alphabet, gap_alphabet)
    {
        return true;
    }
    friend bool operator!=(gap_alphabet, gap_alphabet)
    {
        return false;
    }
};

TEST(debug_stream_test, alphabet)
{
    std::ostringstream o;
    seqan3::debug_stream_type my_stream{o};

    my_stream << gap_alphabet{};

    o.flush();
    EXPECT_EQ(o.str(), "-");
}
