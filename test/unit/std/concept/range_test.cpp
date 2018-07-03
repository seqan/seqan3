// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

#include <gtest/gtest.h>

#include <array>
#include <deque>
#include <forward_list>
#include <list>
#include <string>
#include <vector>

#include <range/v3/view/any_view.hpp>

#include <sdsl/int_vector.hpp>

#include <seqan3/std/concept/range.hpp>

using namespace seqan3;

TEST(range_concepts, range_concept)
{
    EXPECT_TRUE ((range_concept<std::forward_list<char>>));
    EXPECT_TRUE ((range_concept<std::list<char>>));
    EXPECT_TRUE ((range_concept<std::array<char, 2>>));
    EXPECT_TRUE ((range_concept<std::vector<char>>));
    EXPECT_TRUE ((range_concept<std::deque<char>>));
    EXPECT_TRUE ((range_concept<std::string>));

    EXPECT_TRUE ((range_concept<sdsl::int_vector<8>>));
    EXPECT_TRUE ((range_concept<sdsl::int_vector<9>>));

    EXPECT_TRUE ((range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((range_concept<std::vector<char> const>));
    EXPECT_TRUE ((range_concept<std::vector<char> const &>));
}

TEST(range_concepts, input_range_concept)
{
    EXPECT_TRUE ((input_range_concept<std::forward_list<char>>));
    EXPECT_TRUE ((input_range_concept<std::list<char>>));
    EXPECT_TRUE ((input_range_concept<std::array<char, 2>>));
    EXPECT_TRUE ((input_range_concept<std::vector<char>>));
    EXPECT_TRUE ((input_range_concept<std::deque<char>>));
    EXPECT_TRUE ((input_range_concept<std::string>));

    EXPECT_TRUE ((input_range_concept<sdsl::int_vector<8>>));
    EXPECT_TRUE ((input_range_concept<sdsl::int_vector<9>>));

    EXPECT_TRUE ((input_range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((input_range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((input_range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((input_range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((input_range_concept<std::vector<char> const>));
    EXPECT_TRUE ((input_range_concept<std::vector<char> const &>));
}

TEST(range_concepts, forward_range_concept)
{
    EXPECT_TRUE ((forward_range_concept<std::forward_list<char>>));
    EXPECT_TRUE ((forward_range_concept<std::list<char>>));
    EXPECT_TRUE ((forward_range_concept<std::array<char, 2>>));
    EXPECT_TRUE ((forward_range_concept<std::vector<char>>));
    EXPECT_TRUE ((forward_range_concept<std::deque<char>>));
    EXPECT_TRUE ((forward_range_concept<std::string>));

    EXPECT_TRUE ((forward_range_concept<sdsl::int_vector<8>>));
    EXPECT_TRUE ((forward_range_concept<sdsl::int_vector<9>>));

    EXPECT_FALSE((forward_range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((forward_range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((forward_range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((forward_range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((forward_range_concept<std::vector<char> const>));
    EXPECT_TRUE ((forward_range_concept<std::vector<char> const &>));
}

TEST(range_concepts, bidirectional_range_concept)
{
    EXPECT_FALSE((bidirectional_range_concept<std::forward_list<char>>));
    EXPECT_TRUE ((bidirectional_range_concept<std::list<char>>));
    EXPECT_TRUE ((bidirectional_range_concept<std::array<char, 2>>));
    EXPECT_TRUE ((bidirectional_range_concept<std::vector<char>>));
    EXPECT_TRUE ((bidirectional_range_concept<std::deque<char>>));
    EXPECT_TRUE ((bidirectional_range_concept<std::string>));

    EXPECT_TRUE ((bidirectional_range_concept<sdsl::int_vector<8>>));
    EXPECT_TRUE ((bidirectional_range_concept<sdsl::int_vector<9>>));

    EXPECT_FALSE((bidirectional_range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((bidirectional_range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((bidirectional_range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((bidirectional_range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((bidirectional_range_concept<std::vector<char> const>));
    EXPECT_TRUE ((bidirectional_range_concept<std::vector<char> const &>));
}

TEST(range_concepts, random_access_range_concept)
{
    EXPECT_FALSE((random_access_range_concept<std::forward_list<char>>));
    EXPECT_FALSE((random_access_range_concept<std::list<char>>));
    EXPECT_TRUE ((random_access_range_concept<std::array<char, 2>>));
    EXPECT_TRUE ((random_access_range_concept<std::vector<char>>));
    EXPECT_TRUE ((random_access_range_concept<std::deque<char>>));
    EXPECT_TRUE ((random_access_range_concept<std::string>));

    EXPECT_TRUE ((input_range_concept<sdsl::int_vector<8>>));
    EXPECT_TRUE ((input_range_concept<sdsl::int_vector<9>>));

    EXPECT_FALSE((random_access_range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((random_access_range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((random_access_range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((random_access_range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((random_access_range_concept<std::vector<char> const>));
    EXPECT_TRUE ((random_access_range_concept<std::vector<char> const &>));
}

TEST(range_concepts, contiguous_range_concept)
{
    EXPECT_FALSE((contiguous_range_concept<std::forward_list<char>>));
    EXPECT_FALSE((contiguous_range_concept<std::list<char>>));
    EXPECT_TRUE ((contiguous_range_concept<std::array<char, 2>>));
    EXPECT_TRUE ((contiguous_range_concept<std::vector<char>>));
    EXPECT_FALSE((contiguous_range_concept<std::deque<char>>));
    EXPECT_TRUE ((contiguous_range_concept<std::string>));

//     EXPECT_TRUE ((contiguous_range_concept<sdsl::int_vector<8>>));
    EXPECT_FALSE((contiguous_range_concept<sdsl::int_vector<9>>));

    EXPECT_FALSE((contiguous_range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((contiguous_range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((contiguous_range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_FALSE((contiguous_range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((contiguous_range_concept<std::vector<char> const>));
    EXPECT_TRUE ((contiguous_range_concept<std::vector<char> const &>));
}

TEST(range_concepts, output_range_concept)
{
    EXPECT_TRUE ((output_range_concept<std::forward_list<char>, char>));
    EXPECT_TRUE ((output_range_concept<std::list<char>, char>));
    EXPECT_TRUE ((output_range_concept<std::array<char, 2>, char>));
    EXPECT_TRUE ((output_range_concept<std::vector<char>, char>));
    EXPECT_TRUE ((output_range_concept<std::deque<char>, char>));
    EXPECT_TRUE ((output_range_concept<std::string, char>));

    EXPECT_TRUE ((output_range_concept<sdsl::int_vector<8>, uint8_t>));
//     EXPECT_TRUE ((output_range_concept<sdsl::int_vector<9>, uint8_t>));

    EXPECT_FALSE((output_range_concept<ranges::any_view<char, ranges::category::input>, char>));
    EXPECT_FALSE((output_range_concept<ranges::any_view<char, ranges::category::forward>, char>));
    EXPECT_FALSE((output_range_concept<ranges::any_view<char, ranges::category::bidirectional>, char>));
    EXPECT_FALSE((output_range_concept<ranges::any_view<char, ranges::category::random_access>, char>));

    EXPECT_TRUE ((output_range_concept<std::vector<char> &, char>));
    EXPECT_FALSE((output_range_concept<std::vector<char> const, char>));
    EXPECT_FALSE((output_range_concept<std::vector<char> const &, char>));
}

TEST(range_concepts, sized_range_concept)
{
    EXPECT_FALSE((sized_range_concept<std::forward_list<char>>));
    EXPECT_TRUE ((sized_range_concept<std::list<char>>));
    EXPECT_TRUE ((sized_range_concept<std::array<char, 2>>));
    EXPECT_TRUE ((sized_range_concept<std::vector<char>>));
    EXPECT_TRUE ((sized_range_concept<std::deque<char>>));
    EXPECT_TRUE ((sized_range_concept<std::string>));

    EXPECT_TRUE ((sized_range_concept<sdsl::int_vector<8>>));
    EXPECT_TRUE ((sized_range_concept<sdsl::int_vector<9>>));

    EXPECT_FALSE((sized_range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((sized_range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((sized_range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_FALSE((sized_range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((sized_range_concept<std::vector<char> &>));
    EXPECT_TRUE ((sized_range_concept<std::vector<char> const>));
    EXPECT_TRUE ((sized_range_concept<std::vector<char> const &>));
}

TEST(range_concepts, view_concept)
{
    EXPECT_FALSE((view_concept<std::forward_list<char>>));
    EXPECT_FALSE((view_concept<std::list<char>>));
    EXPECT_FALSE((view_concept<std::array<char, 2>>));
    EXPECT_FALSE((view_concept<std::vector<char>>));
    EXPECT_FALSE((view_concept<std::deque<char>>));
    EXPECT_FALSE((view_concept<std::string>));

    EXPECT_FALSE((view_concept<sdsl::int_vector<8>>));
    EXPECT_FALSE((view_concept<sdsl::int_vector<9>>));

    EXPECT_TRUE ((view_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((view_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((view_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((view_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_FALSE((view_concept<std::vector<char> &>));
    EXPECT_FALSE((view_concept<std::vector<char> const>));
    EXPECT_FALSE((view_concept<std::vector<char> const &>));
}

TEST(range_concepts, viewable_range_concept)
{
    EXPECT_FALSE((viewable_range_concept<std::forward_list<char>>));
    EXPECT_FALSE((viewable_range_concept<std::list<char>>));
    EXPECT_FALSE((viewable_range_concept<std::array<char, 2>>));
    EXPECT_FALSE((viewable_range_concept<std::vector<char>>));
    EXPECT_FALSE((viewable_range_concept<std::deque<char>>));
    EXPECT_FALSE((viewable_range_concept<std::string>));

    EXPECT_FALSE((viewable_range_concept<sdsl::int_vector<8>>));
    EXPECT_FALSE((viewable_range_concept<sdsl::int_vector<9>>));

    EXPECT_TRUE ((viewable_range_concept<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((viewable_range_concept<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((viewable_range_concept<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((viewable_range_concept<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((viewable_range_concept<std::vector<char> &>));
    EXPECT_FALSE((viewable_range_concept<std::vector<char> const>));
    EXPECT_TRUE ((viewable_range_concept<std::vector<char> const &>));
}
