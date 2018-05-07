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

#include <sstream>
#include <tuple>
#include <utility>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/quality/all.hpp>

using namespace seqan3;

template <typename T>
class quality : public ::testing::Test
{};

// Q @H-2: testing only {assign,to}_{phred,char,rank} operators, right
// add all alphabets from the quality sub module here
using quality_types     = ::testing::Types<phred42, phred63, phred68>;
using quality_types2    = meta::list<phred42, phred63, phred68>;

TYPED_TEST_CASE(quality, quality_types);

// test assign_char and to_char
TYPED_TEST(quality, assign_to_char)
{
    using char_type = typename TypeParam::char_type;
    // phred42: ['!' .. '^'], phred63: ['!' .. '_'], phred68: [';' .. '~']
    std::vector<char_type> input_char;
    // phred42: ['!' .. 'I', 'J' .. 'J'], phred63: =input_char, phred68: =intput_char
    std::vector<char_type> output_char;

    if constexpr (std::is_same_v<TypeParam, phred42>)
    {
        for (char_type c = '!'; c <= '^'; ++c)
        {
            input_char.push_back(c);
            if (c < 'J')
                output_char.push_back(c);
            else
                output_char.push_back('J');
        }
    }
    else if constexpr (std::is_same_v<TypeParam, phred63>)
    {
        for (char_type c = '!'; c <= '_'; ++c)
            input_char.push_back(c);
        output_char = input_char;
    }
    else if constexpr (std::is_same_v<TypeParam, phred68>)
    {
        for (char_type c = ';'; c <= '~'; ++c)
            input_char.push_back(c);
        output_char = input_char;
    }

    for (auto [ ch, cm ] : ranges::view::zip(input_char, output_char))
        EXPECT_EQ((assign_char(TypeParam{}, ch)).to_char(), cm);
}

// test '=', assign_phred and to_phred
TYPED_TEST(quality, assign_to_phred)
{
    using phred_type = typename TypeParam::phred_type;
    // phred42: ['!'-'!' .. '^'-'!'], phred63: ['!'-'!' .. '_'-'!'], phred68: [';'-';'-5 .. '~'-';'-5]
    std::vector<TypeParam> input_TypeParam_assign_op;
    std::vector<TypeParam> input_TypeParam_assign_fc;
    // phred42: ['!'-'!' .. 'I'-'!', 'J'-'!' .. 'J'-'!'], phred63: =input_phred, phred68: =intput_phred
    std::vector<phred_type> output_phred;

    if constexpr (std::is_same_v<TypeParam, phred42>)
    {
        for (char c = '!'; c <= '^'; ++c)
        {
            TypeParam phred1; phred1 = c - '!';
            input_TypeParam_assign_op.push_back(phred1);
            TypeParam phred2;
            phred2.assign_phred(c - '!');
            input_TypeParam_assign_fc.push_back(phred2);
            if (c < 'J')
                output_phred.push_back(c - '!');
            else
                output_phred.push_back('J' - '!');
        }
    }
    else if constexpr (std::is_same_v<TypeParam, phred63>)
    {
        for (char c = '!'; c <= '_'; ++c){
            TypeParam phred1; phred1 = c - '!';
            input_TypeParam_assign_op.push_back(phred1);
            TypeParam phred2;
            phred2.assign_phred(c - '!');
            input_TypeParam_assign_fc.push_back(phred2);
            output_phred.push_back(c - '!');
        }
    }
    else if constexpr (std::is_same_v<TypeParam, phred68>)
    {
        for (char c = ';'; c <= '~'; ++c){
            TypeParam phred1; phred1 = c - ';' + (-5);
            input_TypeParam_assign_op.push_back(phred1);
            TypeParam phred2;
            phred2.assign_phred(c - ';' + (-5));
            input_TypeParam_assign_fc.push_back(phred2);
            output_phred.push_back(c - ';' + (-5));
        }
    }

    for (auto [ ch, cm ] : ranges::view::zip(input_TypeParam_assign_op, output_phred))
        EXPECT_EQ(ch.to_phred(), cm);

    for (auto [ ch, cm ] : ranges::view::zip(input_TypeParam_assign_fc, output_phred))
        EXPECT_EQ(ch.to_phred(), cm);

}

// test assign_rank and to_rank
TYPED_TEST(quality, assign_to_rank)
{
    using rank_type = typename TypeParam::rank_type; // typename quality_assign_char_Test<gtest_TypeParam_>::TypeParam::char_type;
    // phred42: [0 .. 61], phred63: [0 .. 62], phred68: [0 .. 67]
    std::vector<rank_type> input_rank;
    // phred42: [0 .. 40, 41 .. 41], phred63: =input_rank, phred68: =intput_rank
    std::vector<rank_type> output_rank;

    if constexpr (std::is_same_v<TypeParam, phred42>)
    {
        for (char c = '!'; c <= '^'; ++c)
        {
            input_rank.push_back(c - '!');
            if (c < 'J')
                output_rank.push_back(c - '!');
            else
                output_rank.push_back('J' - '!');
        }
    }
    else if constexpr (std::is_same_v<TypeParam, phred63>)
    {
        for (char c = '!'; c <= '_'; ++c)
            input_rank.push_back(c - '!');
        output_rank = input_rank;
    }
    else if constexpr (std::is_same_v<TypeParam, phred68>)
    {
        for (char c = ';'; c <= '~'; ++c)
            input_rank.push_back(c - ';');
        output_rank = input_rank;
    }

    for (auto [ ch, cm ] : ranges::view::zip(input_rank, output_rank))
        EXPECT_EQ((assign_rank(TypeParam{}, ch)).to_rank(), cm);
}


// ------------------------------------------------------------------
// conversion
// ------------------------------------------------------------------
template<typename Function, typename... Args>
void va_for_each(Function&& fn, Args&&... args)
{
    using dummy = int[];
    static_cast<void>(dummy
        {
            0, (static_cast<void>(fn(std::forward<Args>(args))), 0)...
        });
}

TYPED_TEST(quality, implicit_conversion)
{
    auto action = [](auto other) {
        if constexpr (!std::is_same_v<TypeParam, decltype(other)>)
        {
            TypeParam p;
            p.assign_rank(0);
            [[maybe_unused]] decltype(other) p_other = static_cast<decltype(other)>(p);
            EXPECT_EQ(p_other.to_phred(), decltype(other)::phred_type{0});
        }
    };
    phred42 p42;
    phred63 p63;
    phred68 p68;
    va_for_each(action, p42, p63, p68);
}
