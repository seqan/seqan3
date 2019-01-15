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

TEST(aa27, assign_char)
{
    using t = aa27;
    std::vector<char> chars
    {
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        '*', '!'
    };

    std::vector<aa27> alphabets
    {
        t::A, t::B, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::J, t::K, t::L, t::M,
        t::A, t::B, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::J, t::K, t::L, t::M,
        t::N, t::O, t::P, t::Q, t::R, t::S, t::T, t::U, t::V, t::W, t::X, t::Y, t::Z,
        t::N, t::O, t::P, t::Q, t::R, t::S, t::T, t::U, t::V, t::W, t::X, t::Y, t::Z,
        t::TERMINATOR, t::X
    };

    for (auto [ chr, alp ] : ranges::view::zip(chars, alphabets))
        EXPECT_EQ((assign_char(aa27{}, chr)), alp);
}

TEST(aa27, to_char)
{
    using t = aa27;
    std::vector<char> chars
    {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
        'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'B', 'J', 'O', 'U', 'X', 'Z',
        '*', 'X'
    };

    std::vector<aa27> alphabets
    {
        t::A, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::K, t::L, t::M, t::N, t::P,
        t::Q, t::R, t::S, t::T, t::V, t::W, t::Y, t::B, t::J, t::O, t::U, t::X, t::Z,
        t::TERMINATOR, t::UNKNOWN
    };

    for (auto [ alp, chr ] : ranges::view::zip(alphabets, chars))
        EXPECT_EQ(to_char(alp), chr);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(literals, vector)
{
    aa27_vector v27;
    v27.resize(5, aa27::A);
    EXPECT_EQ(v27, "AAAAA"_aa27);

    std::vector<aa27> w27{aa27::A, aa27::Y, aa27::P, aa27::T, aa27::U, aa27::N, aa27::X, aa27::UNKNOWN,
                        aa27::TERMINATOR};
    EXPECT_EQ(w27, "AYPTUNXX*"_aa27);
}
