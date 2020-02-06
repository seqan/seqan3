// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include "configuration_mock.hpp"

#include <seqan3/core/algorithm/configuration.hpp>

TEST(configuration, concept_check)
{
    EXPECT_TRUE(seqan3::detail::config_element<bar>);
    EXPECT_FALSE(seqan3::detail::config_element<int>);

    EXPECT_TRUE((seqan3::tuple_like<seqan3::configuration<bax, bar>>));
}

TEST(configuration, tuple_size)
{
    EXPECT_EQ((std::tuple_size_v<seqan3::configuration<bax, bar>>), 2u);
    EXPECT_EQ((std::tuple_size<seqan3::configuration<bax, bar>>::value), 2u);
}

TEST(configuration, tuple_element)
{
    EXPECT_TRUE((std::is_same_v<typename std::tuple_element<0, seqan3::configuration<bax, bar>>::type, bax>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, seqan3::configuration<bax, bar>>, bax>));
}

TEST(configuration, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<seqan3::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_copy_constructible_v<seqan3::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_move_constructible_v<seqan3::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_copy_assignable_v<seqan3::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_move_assignable_v<seqan3::configuration<bax, bar>>));
}

TEST(configuration, construction_from_elements)
{
    seqan3::configuration cfg0{};       // empty
    seqan3::configuration cfg1{bax{}};  // one element

    EXPECT_EQ((std::tuple_size_v<decltype(cfg0)>), 0u);
    EXPECT_EQ((std::tuple_size_v<decltype(cfg1)>), 1u);
}

TEST(configuration, size)
{
    seqan3::configuration<foobar<>> cfg{};
    EXPECT_EQ(cfg.size(), 1u);
    EXPECT_EQ((seqan3::configuration<foo, foobar<>>{}.size()), 2u);
    EXPECT_EQ(seqan3::configuration{}.size(), 0u);
}

TEST(configuration, get_by_position)
{
    seqan3::configuration cfg = bax{2.2} | bar{1};

    { // l-value
        EXPECT_EQ(std::get<1>(cfg).value, 1);
        std::get<1>(cfg).value = 3;
        EXPECT_EQ(std::get<1>(cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg)), bar &>));
    }

    { // const l-value
        seqan3::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(std::get<1>(cfg_c).value, 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg_c)), bar const &>));
    }

    { // r-value
        seqan3::configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_r)).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_r))), bar &&>));
    }

    { // const r-value
        seqan3::configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_rc)).value, 3);
        // TODO(rrahn): Enable when seqan3::get(const &&) is fixed for gcc7 as well.
        // EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_rc))), bar const &&>));
    }
}

TEST(configuration, get_by_type)
{
    seqan3::configuration cfg = bax{2.2} | bar{1};

    { // l-value
        EXPECT_FLOAT_EQ(std::get<bax>(cfg).value, 2.2);
        std::get<bax>(cfg).value = 3.1;
        std::get<bar>(cfg).value = 3;
        EXPECT_FLOAT_EQ(std::get<bax>(cfg).value, 3.1);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<bax>(cfg)), bax &>));
    }

    { // const l-value
        seqan3::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(std::get<bar>(cfg_c).value, 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(cfg_c)), bar const &>));
    }

    { // r-value
        seqan3::configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(std::get<bar>(std::move(cfg_r)).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(std::move(cfg_r))), bar &&>));
    }

    { // const r-value
        seqan3::configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(std::get<bar>(std::move(cfg_rc)).value, 3);
        // TODO(rrahn): Enable when seqan3::get(const &&) is fixed for gcc7 as well.
        // EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(std::move(cfg_rc))), bar const &&>));
    }
}

TEST(configuration, get_by_type_template)
{
    seqan3::configuration cfg = bar{1} | foobar<>{std::vector{0, 1, 2, 3}};

    { // l-value
        EXPECT_TRUE(ranges::equal(seqan3::get<foobar>(cfg).value, std::vector{0, 1, 2, 3}));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(cfg)), foobar<> &>));
    }

    { // const l-value
        seqan3::configuration<bar, foobar<>> const cfg_c{cfg};
        EXPECT_TRUE(ranges::equal(seqan3::get<foobar>(cfg_c).value, std::vector{0, 1, 2, 3}));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(cfg_c)), foobar<> const &>));
    }

    { // r-value
        seqan3::configuration<bar, foobar<>> cfg_r{cfg};
        EXPECT_TRUE(ranges::equal(seqan3::get<foobar>(std::move(cfg_r)).value, std::vector{0, 1, 2, 3}));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(std::move(cfg_r))), foobar<> &&>));
    }

    { // const r-value
        seqan3::configuration<bar, foobar<>> const cfg_cr{cfg};
        EXPECT_TRUE(ranges::equal(seqan3::get<foobar>(std::move(cfg_cr)).value, std::vector{0, 1, 2, 3}));
        // TODO(rrahn): Enable when seqan3::get(const &&) is fixed for gcc7 as well.
        // EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(std::move(cfg_cr))), foobar<> const &&>));
    }
}

TEST(configuration, exists_by_type)
{
    seqan3::configuration<bax, bar> cfg{};

    EXPECT_TRUE(std::remove_reference_t<decltype(cfg)>::exists<bax>());
    EXPECT_FALSE(decltype(cfg)::exists<foo>());
}

TEST(configuration, exists_by_type_template)
{
    seqan3::configuration<bax, foobar<>> cfg{};

    EXPECT_TRUE(decltype(cfg)::exists<foobar>());
    EXPECT_TRUE(decltype(cfg)::exists<bax>());
    EXPECT_FALSE(decltype(cfg)::exists<foo>());
}

TEST(configuration, value_or_by_type)
{
    seqan3::configuration cfg = bax{2.2} | bar{1};

    { // l-value
        EXPECT_FLOAT_EQ(cfg.value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(cfg.value_or<foo>(1.3), 1.3);
    }

    { // const l-value
        seqan3::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_FLOAT_EQ(cfg_c.value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(cfg_c.value_or<foo>(1.3), 1.3);
    }

    { // r-value
        seqan3::configuration<bax, bar> cfg_r{cfg};
        EXPECT_FLOAT_EQ(std::move(cfg_r).value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(std::move(cfg_r).value_or<foo>(1.3), 1.3);
    }

    { // const r-value
        seqan3::configuration<bax, bar> const cfg_cr{cfg};
        EXPECT_FLOAT_EQ(std::move(cfg_cr).value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(std::move(cfg_cr).value_or<foo>(1.3), 1.3);
    }
}

TEST(configuration, value_or_by_type_template)
{
    seqan3::configuration cfg = bar{1} | foobar<>{std::vector<int>{0, 1, 2, 3}};

    { // l-value
        EXPECT_TRUE(ranges::equal(cfg.value_or<foobar>(3.3), std::vector{0, 1, 2, 3}));
        EXPECT_FLOAT_EQ(cfg.value_or<foo>(1.3), 1.3);
    }

    { // const l-value
        seqan3::configuration<bar, foobar<>> const cfg_c{cfg};
        EXPECT_TRUE(ranges::equal(cfg_c.value_or<foobar>(3.3), std::vector{0, 1, 2, 3}));
        EXPECT_FLOAT_EQ(cfg_c.value_or<foo>(1.3), 1.3);
    }

    { // r-value
        seqan3::configuration<bar, foobar<>> cfg_r{cfg};
        EXPECT_TRUE(ranges::equal(std::move(cfg_r).value_or<foobar>(3.3), std::vector{0, 1, 2, 3}));
        EXPECT_FLOAT_EQ(std::move(cfg_r).value_or<foo>(1.3), 1.3);
    }

    { // const r-value
        seqan3::configuration<bar, foobar<>> const cfg_cr{cfg};
        EXPECT_TRUE(ranges::equal(std::move(cfg_cr).value_or<foobar>(3.3), std::vector{0, 1, 2, 3}));
        EXPECT_FLOAT_EQ(std::move(cfg_cr).value_or<foo>(1.3), 1.3);
    }
}
