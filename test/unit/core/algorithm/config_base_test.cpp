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

#include <seqan3/core/algorithm/config_base.hpp>
#include <seqan3/core/algorithm/configurator.hpp>

using namespace seqan3;

class bar : public detail::config_base<bar>
{
    friend class detail::config_access<bar>;

    int state{1};
};

TEST(config_base, concept)
{
    EXPECT_TRUE(detail::config_concept<bar>);
    EXPECT_FALSE(detail::config_concept<int>);
}

TEST(config_base, construction)
{
    detail::configurator<bar> cfg{};
    get<0>(cfg) = 3;

    bar b{cfg};
    EXPECT_EQ(b.data(), 3);
}

TEST(config_base, get)
{
    { // l-value
        bar br{};

        EXPECT_EQ(br.data(), 1);
        br.data() = 2;
        EXPECT_EQ(br.data(), 2);

        EXPECT_TRUE((std::is_same_v<decltype(br.data()), int &>));
    }

    { // const l-value
        bar const br_c{};

        EXPECT_EQ(br_c.data(), 1);

        bar br;
        br.data() = 2;
        bar const br_c2{br};

        EXPECT_EQ(br_c2.data(), 2);

        EXPECT_TRUE((std::is_same_v<decltype(br_c2.data()), int const &>));
    }

    { // r-value
        bar br{};

        EXPECT_EQ(std::move(br).data(), 1);
        br.data() = 2;
        EXPECT_EQ(std::move(br).data(), 2);

        EXPECT_TRUE((std::is_same_v<decltype(std::move(br).data()), int &&>));
    }

    { // const r-value
        bar const br_c{};

        EXPECT_EQ(std::move(br_c).data(), 1);

        bar br;
        br.data() = 2;
        bar const br_c2{br};

        EXPECT_EQ(std::move(br_c2).data(), 2);

        EXPECT_TRUE((std::is_same_v<decltype(std::move(br_c2).data()), int const &&>));
    }
}
