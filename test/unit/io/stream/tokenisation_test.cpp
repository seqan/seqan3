// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
// ==========================================================================

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/output_iterator_conversion_adaptor.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/io/stream/tokenisation.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/view/to_char.hpp>

template <typename rng_type>
class tokenisation : public ::testing::Test
{
    virtual void SetUp()
    {
        std::string input = "acgt\tacgt\nacgt acgt\r\nacgtn\n>123;@#\nACGTR\n";
        data = {std::begin(input), std::end(input)};
    }

public:
    rng_type data;
};

// add all <out_rng,in_rng> pairs here.
using input_types = ::testing::Types<std::vector<char>>;

TYPED_TEST_CASE(tokenisation, input_types);

using namespace std::literals;
using namespace seqan3;

TYPED_TEST(tokenisation, transfer_data_w_delim_w_asserter)
{
    std::vector<dna5> target;
    auto target_it = detail::make_conversion_output_iterator(target);
    auto src = this->data | view::single_pass_input;

    detail::transfer_data(target_it, src, is_space, parse_asserter{is_in_alphabet<dna5>{}});
    EXPECT_EQ(std::string{target | view::to_char}, "ACGT"s);

    detail::transfer_data(detail::ignore_output_iterator{}, src, is_alpha, is_space);
    EXPECT_EQ(std::string{target | view::to_char}, "ACGT"s);

    EXPECT_THROW((detail::transfer_data(target_it, src, is_blank, parse_asserter{is_in_alphabet<dna5>{}})), parse_error);
    EXPECT_EQ(std::string{target | view::to_char}, "ACGTACGT"s);
}

TYPED_TEST(tokenisation, transfer_data_w_delim_wo_asserter)
{
    std::vector<dna5> target;
    auto target_it = detail::make_conversion_output_iterator(target);
    auto src = this->data | view::single_pass_input;

    detail::transfer_data(target_it, src, is_space);
    EXPECT_EQ(std::string{target | view::to_char}, "ACGT"s);

    detail::transfer_data(detail::ignore_output_iterator{}, src, is_alpha);
    EXPECT_EQ(std::string{target | view::to_char}, "ACGT"s);

    EXPECT_NO_THROW((detail::transfer_data(target_it, src, is_blank)));
    EXPECT_EQ(std::string{target | view::to_char}, "ACGTACGTNACGT"s);
}

TYPED_TEST(tokenisation, read_until)
{
    std::vector<dna5> target;
    auto target_it = detail::make_conversion_output_iterator(target);
    auto src = this->data | view::single_pass_input;
    parse_asserter asserter{is_in_alphabet<dna5>{}};

    read_until(target_it, src, is_space, asserter);
    EXPECT_EQ(std::string{target | view::to_char}, "ACGT"s);

    unsigned count = 1;
    read_until(std::ignore, src, [&](auto){ return (count-- > 0) ? false : true; });
    EXPECT_EQ(std::string{target | view::to_char}, "ACGT"s);

    EXPECT_THROW((read_until(target_it, src, is_blank, asserter)), parse_error);
    EXPECT_EQ(std::string{target | view::to_char}, "ACGTACGT"s);
}
