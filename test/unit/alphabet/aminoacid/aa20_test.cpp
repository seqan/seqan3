// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <range/v3/view/zip.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"
#include "aminoacid_test_template.hpp"

#include <seqan3/alphabet/aminoacid/aa20.hpp>

INSTANTIATE_TYPED_TEST_CASE_P(aa20, alphabet, aa20);
INSTANTIATE_TYPED_TEST_CASE_P(aa20, alphabet_constexpr, aa20);
INSTANTIATE_TYPED_TEST_CASE_P(aa20, aminoacid, aa20);

TEST(aa20, assign_char)
{
    using t = aa20;
    std::vector<char> chars
    {
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        '*', '!'
    };

    std::vector<aa20> alphabets
    {
        t::A, t::D, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::L, t::K, t::L, t::M,
        t::A, t::D, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::L, t::K, t::L, t::M,
        t::N, t::L, t::P, t::Q, t::R, t::S, t::T, t::C, t::V, t::W, t::S, t::Y, t::E,
        t::N, t::L, t::P, t::Q, t::R, t::S, t::T, t::C, t::V, t::W, t::S, t::Y, t::E,
        t::W, t::S
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ((assign_char(aa20{}, chr)), alp);
}

TEST(aa20, to_char)
{
    using t = aa20;
    std::vector<char> chars
    {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
        'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'D', 'L', 'L', 'C', 'S', 'E',
        'W', 'S'
    };

    std::vector<aa20> alphabets
    {
        t::A, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::K, t::L, t::M, t::N, t::P,
        t::Q, t::R, t::S, t::T, t::V, t::W, t::Y, t::B, t::J, t::O, t::U, t::X, t::Z,
        t::TERMINATOR, t::UNKNOWN
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ(to_char(alp), chr);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(literals, vector)
{
    aa20_vector v20;
    v20.resize(5, aa20::B);
    EXPECT_EQ(v20, "DDDDD"_aa20);

    std::vector<aa20> w20{aa20::A, aa20::B, aa20::J, aa20::O, aa20::U, aa20::X, aa20::Z, aa20::UNKNOWN,
                        aa20::TERMINATOR, aa20::TERMINATOR};
    EXPECT_EQ(w20, "ADLLCSESW*"_aa20);
}
