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

#include <seqan3/core/algorithm/configuration.hpp>

using ::testing::ElementsAre;

using namespace seqan3;

TEST(configuration, concept_check)
{
    EXPECT_TRUE(detail::config_element_concept<bar>);
    EXPECT_FALSE(detail::config_element_concept<int>);

    EXPECT_TRUE((tuple_like_concept<configuration<bax, bar>>));
}

TEST(configuration, tuple_size)
{
    EXPECT_EQ((std::tuple_size_v<configuration<bax, bar>>), 2u);
    EXPECT_EQ((std::tuple_size<configuration<bax, bar>>::value), 2u);
}

TEST(configuration, tuple_element)
{
    EXPECT_TRUE((std::is_same_v<typename std::tuple_element<0, configuration<bax, bar>>::type, bax>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, configuration<bax, bar>>, bax>));
}

TEST(configuration, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<configuration<bax, bar>>));
    EXPECT_TRUE((std::is_copy_constructible_v<configuration<bax, bar>>));
    EXPECT_TRUE((std::is_move_constructible_v<configuration<bax, bar>>));
    EXPECT_TRUE((std::is_copy_assignable_v<configuration<bax, bar>>));
    EXPECT_TRUE((std::is_move_assignable_v<configuration<bax, bar>>));
}

TEST(configuration, construction_from_elements)
{
    configuration cfg0{};         // empty
    configuration cfg1{bax{}};      // one element

    EXPECT_EQ((std::tuple_size_v<decltype(cfg0)>), 0u);
    EXPECT_EQ((std::tuple_size_v<decltype(cfg1)>), 1u);
}

TEST(configuration, size)
{
    configuration<foobar<>> cfg{};
    EXPECT_EQ(cfg.size(), 1u);
    EXPECT_EQ((configuration<foo, foobar<>>{}.size()), 2u);
    EXPECT_EQ(configuration{}.size(), 0u);
}

TEST(configuration, get_by_position)
{
    configuration cfg = bax{2.2} | bar{1};

    { // l-value
        EXPECT_EQ(get<1>(cfg).value, 1);
        get<1>(cfg).value = 3;
        EXPECT_EQ(get<1>(cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<1>(cfg)), bar &>));
    }

    { // const l-value
        configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(get<1>(cfg_c).value, 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<1>(cfg_c)), bar const &>));
    }

    { // r-value
        configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(get<1>(std::move(cfg_r)).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<1>(std::move(cfg_r))), bar &&>));
    }

    { // const r-value
        configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(get<1>(std::move(cfg_rc)).value, 3);
        // TODO(rrahn): Enable when get(const &&) is fixed for gcc7 as well.
        // EXPECT_TRUE((std::is_same_v<decltype(get<1>(std::move(cfg_rc))), bar const &&>));
    }
}

TEST(configuration, get_by_type)
{
    configuration cfg = bax{2.2} | bar{1};

    { // l-value
        EXPECT_FLOAT_EQ(get<bax>(cfg).value, 2.2);
        get<bax>(cfg).value = 3.1;
        get<bar>(cfg).value = 3;
        EXPECT_FLOAT_EQ(get<bax>(cfg).value, 3.1);
        EXPECT_TRUE((std::is_same_v<decltype(get<bax>(cfg)), bax &>));
    }

    { // const l-value
        configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(get<bar>(cfg_c).value, 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(cfg_c)), bar const &>));
    }

    { // r-value
        configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(get<bar>(std::move(cfg_r)).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(std::move(cfg_r))), bar &&>));
    }

    { // const r-value
        configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(get<bar>(std::move(cfg_rc)).value, 3);
        // TODO(rrahn): Enable when get(const &&) is fixed for gcc7 as well.
        // EXPECT_TRUE((std::is_same_v<decltype(get<bar>(std::move(cfg_rc))), bar const &&>));
    }
}

TEST(configuration, get_by_type_template)
{
    configuration cfg = bar{1} | foobar<>{std::vector{0, 1, 2, 3}};

    { // l-value
        EXPECT_THAT(get<foobar>(cfg).value, ElementsAre(0, 1, 2, 3));
        EXPECT_TRUE((std::is_same_v<decltype(get<foobar>(cfg)), foobar<> &>));
    }

    { // const l-value
        configuration<bar, foobar<>> const cfg_c{cfg};
        EXPECT_THAT(get<foobar>(cfg_c).value, ElementsAre(0, 1, 2, 3));
        EXPECT_TRUE((std::is_same_v<decltype(get<foobar>(cfg_c)), foobar<> const &>));
    }

    { // r-value
        configuration<bar, foobar<>> cfg_r{cfg};
        EXPECT_THAT(get<foobar>(std::move(cfg_r)).value, ElementsAre(0, 1, 2, 3));
        EXPECT_TRUE((std::is_same_v<decltype(get<foobar>(std::move(cfg_r))), foobar<> &&>));
    }

    { // const r-value
        configuration<bar, foobar<>> const cfg_cr{cfg};
        EXPECT_THAT(get<foobar>(std::move(cfg_cr)).value, ElementsAre(0, 1, 2, 3));
        // TODO(rrahn): Enable when get(const &&) is fixed for gcc7 as well.
        // EXPECT_TRUE((std::is_same_v<decltype(get<foobar>(std::move(cfg_cr))), foobar<> const &&>));
    }
}

TEST(configuration, exists_by_type)
{
    configuration<bax, bar> cfg{};

    EXPECT_TRUE(std::remove_reference_t<decltype(cfg)>::exists<bax>());
    EXPECT_FALSE(decltype(cfg)::exists<foo>());
}

TEST(configuration, exists_by_type_template)
{
    configuration<bax, foobar<>> cfg{};

    EXPECT_TRUE(decltype(cfg)::exists<foobar>());
    EXPECT_TRUE(decltype(cfg)::exists<bax>());
    EXPECT_FALSE(decltype(cfg)::exists<foo>());
}

TEST(configuration, value_or_by_type)
{
    configuration cfg = bax{2.2} | bar{1};

    { // l-value
        EXPECT_FLOAT_EQ(cfg.value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(cfg.value_or<foo>(1.3), 1.3);
    }

    { // const l-value
        configuration<bax, bar> const cfg_c{cfg};
        EXPECT_FLOAT_EQ(cfg_c.value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(cfg_c.value_or<foo>(1.3), 1.3);
    }

    { // r-value
        configuration<bax, bar> cfg_r{cfg};
        EXPECT_FLOAT_EQ(std::move(cfg_r).value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(std::move(cfg_r).value_or<foo>(1.3), 1.3);
    }

    { // const r-value
        configuration<bax, bar> const cfg_cr{cfg};
        EXPECT_FLOAT_EQ(std::move(cfg_cr).value_or<bax>(1.3), 2.2);
        EXPECT_FLOAT_EQ(std::move(cfg_cr).value_or<foo>(1.3), 1.3);
    }
}

TEST(configuration, value_or_by_type_template)
{
    configuration cfg = bar{1} | foobar<>{std::vector<int>{0, 1, 2, 3}};

    { // l-value
        EXPECT_THAT(cfg.value_or<foobar>(3.3), ElementsAre(0, 1, 2, 3));
        EXPECT_FLOAT_EQ(cfg.value_or<foo>(1.3), 1.3);
    }

    { // const l-value
        configuration<bar, foobar<>> const cfg_c{cfg};
        EXPECT_THAT(cfg_c.value_or<foobar>(3.3), ElementsAre(0, 1, 2, 3));
        EXPECT_FLOAT_EQ(cfg_c.value_or<foo>(1.3), 1.3);
    }

    { // r-value
        configuration<bar, foobar<>> cfg_r{cfg};
        EXPECT_THAT(std::move(cfg_r).value_or<foobar>(3.3), ElementsAre(0, 1, 2, 3));
        EXPECT_FLOAT_EQ(std::move(cfg_r).value_or<foo>(1.3), 1.3);
    }

    { // const r-value
        configuration<bar, foobar<>> const cfg_cr{cfg};
        EXPECT_THAT(std::move(cfg_cr).value_or<foobar>(3.3), ElementsAre(0, 1, 2, 3));
        EXPECT_FLOAT_EQ(std::move(cfg_cr).value_or<foo>(1.3), 1.3);
    }
}
