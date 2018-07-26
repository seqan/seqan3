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

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>

using namespace seqan3;

TEST(template_inspect, transfer_template_args_onto_t)
{
    using tl = type_list<int, char, double>;
    using t = detail::transfer_template_args_onto<tl, std::tuple>::type;
    EXPECT_TRUE((std::is_same_v<t, std::tuple<int, char, double>>));

    // shortcut
    using t = detail::transfer_template_args_onto_t<tl, std::tuple>;
    EXPECT_TRUE((std::is_same_v<t, std::tuple<int, char, double>>));
}

TEST(template_inspect, is_type_specialisation_of_v)
{
    using tl = type_list<int, char, double>;
    EXPECT_TRUE((detail::is_type_specialisation_of_v<tl, type_list>));
    EXPECT_FALSE((detail::is_type_specialisation_of_v<int, type_list>));
}

template <int i, char c>
struct t1 {};

template <int _i, char _c>
struct t2
{
    static constexpr auto i = _i;
    static constexpr auto c = _c;
};

TEST(template_inspect, transfer_template_vargs_onto_t)
{
    using tl = t1<1, 'a'>;
    using ta = detail::transfer_template_vargs_onto<tl, t2>::type;
    EXPECT_EQ(1,   ta::i);
    EXPECT_EQ('a', ta::c);

    // shortcut
    using ta2 = detail::transfer_template_vargs_onto_t<tl, t2>;
    EXPECT_EQ(1,   ta2::i);
    EXPECT_EQ('a', ta2::c);
}

TEST(template_inspect, is_value_specialisation_of_v)
{
    using tl = t1<1, 'a'>;

    EXPECT_TRUE((detail::is_value_specialisation_of_v<tl, t1>));
    EXPECT_FALSE((detail::is_value_specialisation_of_v<int, t1>));
}
