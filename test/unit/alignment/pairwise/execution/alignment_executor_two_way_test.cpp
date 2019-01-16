// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
