// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "configuration_mock.hpp"

#include <seqan3/core/algorithm/configuration.hpp>
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

TEST(pipeable_config_element, configuration_with_element)
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

TEST(pipeable_config_element, configuration_with_configuration)
{
    configuration<bar> tmp{};
    configuration<bax> b2{};
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
        auto cfg = tmp | configuration<bax>{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), configuration<bar, bax>>));
    }

    // rvalue | rvalue
    {
        auto cfg = configuration<bar>{} | configuration<bax>{};
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
