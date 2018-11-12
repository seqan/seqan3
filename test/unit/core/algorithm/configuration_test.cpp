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

#include <seqan3/core/algorithm/configuration.hpp>

using namespace seqan3;

struct bar
{
    int value{1};
};

struct bax
{
    float value{2.2};
};

TEST(configuration, concept_check)
{
    EXPECT_TRUE(detail::ConfigElement<bar>);
    EXPECT_FALSE(detail::ConfigElement<int>);

    EXPECT_TRUE((TupleLike<detail::configuration<bax, bar>>));
}

TEST(configuration, tuple_size)
{
    EXPECT_EQ((std::tuple_size_v<detail::configuration<bax, bar>>), 2);
    EXPECT_EQ((std::tuple_size<detail::configuration<bax, bar>>::value), 2);
}

TEST(configuration, tuple_element)
{
    EXPECT_TRUE((std::is_same_v<typename std::tuple_element<0, detail::configuration<bax, bar>>::type, bax>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, detail::configuration<bax, bar>>, bax>));
}

TEST(configuration, metafunction)
{
    EXPECT_TRUE((detail::is_algorithm_configuration<detail::configuration<bax, bar>>::value));
    EXPECT_TRUE((detail::is_algorithm_configuration_v<detail::configuration<bax, bar>>));
    EXPECT_FALSE((detail::is_algorithm_configuration<type_list<bax>>::value));
    EXPECT_FALSE((detail::is_algorithm_configuration_v<type_list<bax>>));
}

TEST(configuration, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<detail::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_copy_constructible_v<detail::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_move_constructible_v<detail::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_copy_assignable_v<detail::configuration<bax, bar>>));
    EXPECT_TRUE((std::is_move_assignable_v<detail::configuration<bax, bar>>));
}

TEST(configuration, get_by_position)
{
    detail::configuration<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(std::get<1>(cfg).value, 1);
        std::get<1>(cfg).value = 3;
        EXPECT_EQ(std::get<1>(cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg)), bar &>));
    }

    { // const l-value
        detail::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(std::get<1>(cfg_c).value, 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg_c)), bar const &>));
    }

    { // r-value
        detail::configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_r)).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_r))), bar &&>));
    }

    { // const r-value
        detail::configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_rc)).value, 3);
        // TODO(rrahn): Enable when std::get(const &&) is fixed for gcc7 as well.
        //EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_rc))), bar const &&>));
    }
}

TEST(configuration, construction_from_tuple)
{
    detail::configuration cfg{std::tuple<bar, bax>{bar{}, bax{}}};

    using t = typename decltype(cfg)::base_type;
    EXPECT_EQ(std::tuple_size_v<decltype(static_cast<t>(cfg))>, 2);
}

template <size_t I>
struct foo
{
    size_t value{I};
};

TEST(configuration, push_front)
{
    { // l-value
        detail::configuration<foo<0>> cfg{};
        auto cfg_2 = cfg.push_front(foo<1>{});

        EXPECT_TRUE((std::is_same_v<decltype(cfg_2), detail::configuration<foo<1>, foo<0>>>));
    }

    { // r-value
        auto cfg_2 = detail::configuration<foo<0>>{}.push_front(foo<1>{});

        EXPECT_TRUE((std::is_same_v<decltype(cfg_2), detail::configuration<foo<1>, foo<0>>>));
    }
}

TEST(configuration, replace_with)
{
    { // l-value
        detail::configuration<foo<0>> cfg{};
        auto cfg_2 = cfg.replace_with(foo<0>{}, foo<1>{});
        EXPECT_TRUE((std::is_same_v<decltype(cfg_2), detail::configuration<foo<1>>>));
    }

    { // r-value
        auto cfg_2 = detail::configuration<foo<0>>{}.replace_with(foo<0>{}, foo<1>{});
        EXPECT_TRUE((std::is_same_v<decltype(cfg_2), detail::configuration<foo<1>>>));
    }
}

TEST(configuration, size)
{
    detail::configuration<foo<0>> cfg{};
    EXPECT_EQ(cfg.size(), 1);
    EXPECT_EQ((detail::configuration<foo<1>, foo<0>>{}.size()), 2);
    EXPECT_EQ(detail::configuration<>{}.size(), 0);
}

struct bar_fn_impl : public detail::configuration_fn_base<bar_fn_impl>
{

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg,
                                           int new_v) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(bar{new_v});
    }

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(bar{0});
    }
};

struct bax_fn_impl : public detail::configuration_fn_base<bax_fn_impl>
{
    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg,
                          float new_v) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(bax{new_v});
    }

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(bax{0.});
    }
};

TEST(configuration, template_deduction_from_proxy)
{
    detail::configuration cfg = bar_fn_impl{}(3);

    EXPECT_EQ(std::get<0>(cfg).value, 3);
    EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar>>));
}

TEST(configuration, template_deduction_from_variable)
{
    detail::configuration cfg = bar_fn_impl{};

    EXPECT_EQ(std::get<0>(cfg).value, 0);
    EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar>>));
}

TEST(configuration_fn, invoke_w_configuration)
{

    { // as l-value
        detail::configuration<bax> cfg;
        auto new_cfg = bar_fn_impl{}(cfg, 3);

        EXPECT_EQ(std::get<0>(new_cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(new_cfg), detail::configuration<bar, bax>>));
    }

    { // as r-value
        auto new_cfg = bar_fn_impl{}(detail::configuration<bax>{}, 3);

        EXPECT_EQ(std::get<0>(new_cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(new_cfg), detail::configuration<bar, bax>>));
    }
}

TEST(configuration_fn, pipeable_w_derived_fn)
{
    { // as l-value
        detail::configuration<bax> cfg{};
        bar_fn_impl fn{};
        auto cfg_r = cfg | fn;

        EXPECT_EQ(std::get<0>(cfg_r).value, 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configuration<bar, bax>>));
    }

    { // as l-value
        bar_fn_impl fn{};
        auto cfg_r = detail::configuration<bax>{} | fn;

        EXPECT_EQ(std::get<0>(cfg_r).value, 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configuration<bar, bax>>));
    }

    { // as r-value
        detail::configuration<bax> cfg{};
        auto cfg_r = cfg | bar_fn_impl{};

        EXPECT_EQ(std::get<0>(cfg_r).value, 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configuration<bar, bax>>));
    }

    { // as r-value
        auto cfg_r = detail::configuration<bax>{} | bar_fn_impl{};

        EXPECT_EQ(std::get<0>(cfg_r).value, 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configuration<bar, bax>>));
    }
}

TEST(configuration_fn, pipeable_w_proxy)
{
    { // as l-value
        detail::configuration<bax> cfg_{};
        int val = 3;
        auto proxy = bar_fn_impl{}(val);
        auto cfg = cfg_ | proxy;

        EXPECT_EQ(std::get<0>(cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }

    { // as l-value
        int val = 3;
        auto proxy = bar_fn_impl{}(val);
        auto cfg = detail::configuration<bax>{} | proxy;

        EXPECT_EQ(std::get<0>(cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }

    { // as r-value
        detail::configuration<bax> cfg_{};
        auto cfg = cfg_ | bar_fn_impl{}(3);

        EXPECT_EQ(std::get<0>(cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }

    { // as r-value
        auto cfg = detail::configuration<bax>{} | bar_fn_impl{}(3);

        EXPECT_EQ(std::get<0>(cfg).value, 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }
}

TEST(configuration_fn, pipable_fn_fn)
{
    auto cfg = bar_fn_impl{} | bax_fn_impl{};

    EXPECT_EQ(std::get<1>(cfg).value, 0);
    EXPECT_FLOAT_EQ(std::get<0>(cfg).value, 0.0);
}

TEST(configuration_fn, pipable_fn_proxy)
{
    auto cfg = bar_fn_impl{} | bax_fn_impl{}(3.0);

    EXPECT_EQ(std::get<1>(cfg).value, 0);
    EXPECT_FLOAT_EQ(std::get<0>(cfg).value, 3.0);
}

TEST(configuration_fn, pipable_proxy_fn)
{
    auto cfg = bar_fn_impl{}(2) | bax_fn_impl{};

    EXPECT_EQ(std::get<1>(cfg).value, 2);
    EXPECT_FLOAT_EQ(std::get<0>(cfg).value, 0.0);
}

TEST(configuration_fn, pipable_proxy_proxy)
{
    auto cfg = bar_fn_impl{}(2) | bax_fn_impl{}(3.0);

    EXPECT_EQ(std::get<1>(cfg).value, 2);
    EXPECT_FLOAT_EQ(std::get<0>(cfg).value, 3.0);
}

TEST(configuration_fn, pipeable_w_empty_config)
{
    {
        detail::configuration cfg;
        auto cfg2 = cfg | bar_fn_impl{}(2);
        EXPECT_EQ(std::get<0>(cfg2).value, 2);
    }

    {
        detail::configuration const cfg;
        auto cfg2 = cfg | bar_fn_impl{}(2);
        EXPECT_EQ(std::get<0>(cfg2).value, 2);
    }
}
