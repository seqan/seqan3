// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
 * \ingroup alphabet
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Tests for the structure alphabets.
 */

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/structure/dot_bracket3.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class structure : public ::testing::Test
{};

// add all alphabets from the structure sub module here
using structure_types  = ::testing::Types<db3>;

TYPED_TEST_CASE(structure, structure_types);

TYPED_TEST(structure, assign_char)
{
    using t = TypeParam;
    std::vector<char> input
    {
        '.', '(', ')',
        'a', 'c', 'g', 't', 'u', 'n'
    };

    std::vector<TypeParam> cmp;
    if constexpr (std::is_same_v<TypeParam, db3>)
    {
        cmp =
        {
        t::NP, t::BL, t::BR,
        t::NA, t::NA, t::NA, t::NA, t::NA, t::NA
        };
    }

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
    EXPECT_EQ((assign_char(TypeParam{}, ch)), cm);
}

TYPED_TEST(structure, to_char)
{
    if constexpr (std::is_same_v<TypeParam, db3> /* || std::is_same_v<TypeParam, wuss18> */ )
    {
        EXPECT_EQ(to_char(TypeParam::NP), '.');
        EXPECT_EQ(to_char(TypeParam::BL), '(');
        EXPECT_EQ(to_char(TypeParam::BR), ')');
        EXPECT_EQ(to_char(TypeParam::NA), '.');
    }
}

TYPED_TEST(structure, stream_operator)
{
    std::stringstream ss;
    ss << TypeParam::BL << TypeParam::BR << TypeParam::NP;
    EXPECT_EQ(ss.str(), "().");
}

TYPED_TEST(structure, concept)
{
    EXPECT_TRUE(structure_concept<TypeParam>);
}

// ------------------------------------------------------------------
// conversion
// ------------------------------------------------------------------
/*
// conversion db3 - wuss18
TYPED_TEST(structure, implicit_conversion)
{
    using complement_type = std::conditional_t<std::is_same_v<TypeParam, db3>, wuss18,
                            std::conditional_t<std::is_same_v<TypeParam, wuss18>, db3,
                            void>>;
    if constexpr (!std::is_same_v<complement_type, void>)
    {
        // construct
        EXPECT_EQ(complement_type{TypeParam::BL}, complement_type::BL);
        // assign
        complement_type l{};
        l = TypeParam::BL;
        EXPECT_EQ(l, complement_type::BL);
    }
}
*/

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(db3_literals, vector)
{
    db3_vector v;
    v.resize(5, db3::BL);
    EXPECT_EQ(v, "((((("_db3);

    std::vector<db3> w{db3::NP, db3::BL, db3::BL, db3::BR, db3::BR, db3::NA};
    EXPECT_EQ(w, ".(())."_db3);
}

TEST(db3_literals, basic_string)
{
    db3_string v;
    v.resize(5, db3::BR);
    EXPECT_EQ(v, ")))))"_db3s);

    std::basic_string<db3, std::char_traits<db3>> w{db3::BL, db3::NP, db3::BR, db3::BL, db3::BR, db3::NA};
    EXPECT_EQ(w, "(.)()."_db3s);
}
