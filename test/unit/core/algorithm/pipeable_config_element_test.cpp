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

#include "configuration_mock.hpp"

#include <seqan3/core/algorithm/pipeable_config_element.hpp>

using namespace seqan3;

TEST(pipeable_config_element, two_elements)
{
    // lvalue | lvalue
    bar b1{};
    bax b2{};
    {
        auto cfg = b1 | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }

    // rvalue | lvalue
    {
        auto cfg = bar{} | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }

    // lvalue | rvalue
    {
        auto cfg = b1 | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }

    // rvalue | rvalue
    {
        auto cfg = bar{} | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }
}

TEST(pipeable_config_element, configuration_with_configuration)
{
    configuration<bar> tmp{};
    bax b2{};
    // lvalue | lvalue
    {
        auto cfg = tmp | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }

    // rvalue | lvalue
    {
        auto cfg = configuration<bar>{} | b2;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }

    // lvalue | rvalue
    {
        auto cfg = tmp | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }

    // rvalue | rvalue
    {
        auto cfg = configuration<bar>{} | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }
}

TEST(pipeable_config_element, multiple_elements)
{
    configuration<bar> tmp{};
    bax b2{};
    {
        auto cfg = foo{} | bar{} | bax{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<foo, bar, bax>>));
    }

    {
        auto cfg = configuration<bar>{} | b2 | foo{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax, foo>>));
    }
}

TEST(pipeable_config_element, const_config)
{
    configuration<foobar<>> const tmp{};

    {
        auto cfg = tmp | foo{} | bar{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<foobar<>, foo, bar>>));
    }
}
