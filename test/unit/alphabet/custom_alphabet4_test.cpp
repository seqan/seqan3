// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>

struct non_writeable_alphabet
{
    using rank_t = uint8_t;
    using char_t = char;

    rank_t to_rank() const noexcept;
    char_t to_char() const noexcept;

    static constexpr bool alphabet_size{1};

    friend bool operator<(non_writeable_alphabet, non_writeable_alphabet);
    friend bool operator<=(non_writeable_alphabet, non_writeable_alphabet);
    friend bool operator>(non_writeable_alphabet, non_writeable_alphabet);
    friend bool operator>=(non_writeable_alphabet, non_writeable_alphabet);
    friend bool operator==(non_writeable_alphabet, non_writeable_alphabet);
    friend bool operator!=(non_writeable_alphabet, non_writeable_alphabet);
};

// see issue https://github.com/seqan/seqan3/issues/1518
TEST(non_writeable_alphabet_test, issue1518)
{
    EXPECT_TRUE(seqan3::semialphabet<non_writeable_alphabet>);
    EXPECT_TRUE(seqan3::alphabet<non_writeable_alphabet>);

    EXPECT_FALSE(seqan3::writable_semialphabet<non_writeable_alphabet>);
    EXPECT_FALSE(seqan3::writable_alphabet<non_writeable_alphabet>);
}
