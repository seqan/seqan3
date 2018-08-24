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

#include <seqan3/io/detail/output_iterator_conversion_adaptor.hpp>
#include <seqan3/range/view/to_char.hpp>

#include <iostream>
#include <utility>
#include <sstream>
#include <string>
#include <vector>

#include <range/v3/utility/iterator.hpp>

using namespace seqan3;
using namespace std::literals;

template <typename T>
class oiter_conversion_adaptor : public ::testing::Test
{
public:

    using output_iterator = std::tuple_element_t<0, T>;
    using output_value_type = std::tuple_element_t<2, T>;
    using output_iterator_adaptor = detail::output_iterator_conversion_adaptor<output_iterator, output_value_type>;

    auto get_output()
    {
        if constexpr (std::is_base_of_v<std::ios_base, std::tuple_element_t<1, T>>)
        {
            return out.str();
        }
        else
        {
            auto v = out | view::to_char;
            return std::string{v.begin(), v.end()};
        }
    }

    std::tuple_element_t<1, T> out;
    output_iterator inner_it{out};
};

using output_value_types = ::testing::Types<std::tuple<ranges::ostream_iterator<char>, std::ostringstream, char>,
                                            std::tuple<ranges::ostreambuf_iterator<char>, std::ostringstream, char>,
                                            std::tuple<ranges::back_insert_iterator<std::vector<char>>,
                                                       std::vector<char>,
                                                       char>,
                                            std::tuple<ranges::back_insert_iterator<std::vector<dna4>>,
                                                       std::vector<dna4>,
                                                       dna4>
                                           >;

TYPED_TEST_CASE(oiter_conversion_adaptor, output_value_types);

TYPED_TEST(oiter_conversion_adaptor, concept)
{
    using adaptor_type = typename TestFixture::output_iterator_adaptor;

    EXPECT_TRUE((std::OutputIterator<adaptor_type, char>));
    EXPECT_FALSE((std::InputIterator<adaptor_type>));
}

TYPED_TEST(oiter_conversion_adaptor, assign)
{
    using adaptor_type = typename TestFixture::output_iterator_adaptor;

    adaptor_type it{this->inner_it};

    *it = 'A';
    EXPECT_EQ(this->get_output(), "A"s);
    it = 'C';
    EXPECT_EQ(this->get_output(), "AC"s);
}

TYPED_TEST(oiter_conversion_adaptor, pre_increment)
{
    using adaptor_type = typename TestFixture::output_iterator_adaptor;

    [[maybe_unused]] adaptor_type it{this->inner_it};

    EXPECT_TRUE((std::Same<std::remove_reference_t<decltype(++it)>, adaptor_type>));
}

TYPED_TEST(oiter_conversion_adaptor, post_increment)
{
    using adaptor_type = typename TestFixture::output_iterator_adaptor;

    [[maybe_unused]] adaptor_type it{this->inner_it};

    EXPECT_TRUE((std::Same<std::remove_reference_t<decltype(it++)>, adaptor_type>));
}

TYPED_TEST(oiter_conversion_adaptor, dereference)
{
    using adaptor_type = typename TestFixture::output_iterator_adaptor;

    [[maybe_unused]] adaptor_type it{this->inner_it};

    EXPECT_TRUE((std::Same<std::remove_reference_t<decltype(*it)>, adaptor_type>));
}

TEST(output_iterator, vector)
{
    std::string input{"12345 6789"};
    std::vector<char> vec;
    auto it = detail::make_conversion_output_iterator(vec);
    for (auto val : input)
    {
        *it = val;
        ++it;
    }
    EXPECT_EQ(vec[0], '1');
    EXPECT_EQ(vec[1], '2');
    EXPECT_EQ(vec[2], '3');
    EXPECT_EQ(vec[3], '4');
    EXPECT_EQ(vec[4], '5');
    EXPECT_EQ(vec[5], ' ');
    EXPECT_EQ(vec[6], '6');
    EXPECT_EQ(vec[7], '7');
    EXPECT_EQ(vec[8], '8');
    EXPECT_EQ(vec[9], '9');
}

TEST(output_iterator, dna_vector)
{
    std::string input{"ACGT TGCA"};
    std::vector<dna4> vec;
    auto it = detail::make_conversion_output_iterator(vec);
    for (auto val : input)
    {
        *it = val;
        ++it;
    }
    EXPECT_EQ(vec[0], dna4{dna4::A});
    EXPECT_EQ(vec[1], dna4{dna4::C});
    EXPECT_EQ(vec[2], dna4{dna4::G});
    EXPECT_EQ(vec[3], dna4{dna4::T});
    EXPECT_EQ(vec[4], dna4{dna4::A});
    EXPECT_EQ(vec[5], dna4{dna4::T});
    EXPECT_EQ(vec[6], dna4{dna4::G});
    EXPECT_EQ(vec[7], dna4{dna4::C});
    EXPECT_EQ(vec[8], dna4{dna4::A});
}

TEST(output_iterator, ostream)
{
    std::string input{"12345 6789"};
    std::ostringstream stream;
    auto it = detail::make_conversion_output_iterator(stream);
    for (auto val : input)
    {
        *it = val;
        ++it;
    }

    EXPECT_EQ(stream.str(), "12345 6789");
}
