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

#include <string>

#include <range/v3/view/bounded.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alignment/pairwise/execution/alignment_executor_two_way.hpp>

using namespace seqan3;

struct foo
{
    int i;

    std::string operator()(std::string & /*ignore*/) const
    {
        return std::string{"my_test_"} + std::to_string(i);
    }
};

struct foo_selector
{
    using result_type = std::string;

    template <typename task_t>
    auto select(task_t && t)
    {
        std::function<result_type(result_type &)> f = t;
        return f;
    }
};

//TODO: Currently we only do a generic integration test. We need to add more testing in the end.
TEST(alignment_excecutor_two_way, stream_integration)
{
    std::vector<foo> resource_rng{foo{0}, foo{1}, foo{2}, foo{3}, foo{4}, foo{0}, foo{1}, foo{2}, foo{3}, foo{4}};

    detail::alignment_executor_two_way exec{resource_rng, foo_selector{}};

    int c = 0;
    for (auto res : exec.range())
    {
        std::string test = "my_test_" + std::to_string(c % 5);
        EXPECT_EQ(test, res);
        ++c;
    }
}
