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

#include <seqan3/core/algorithm/configurator.hpp>
#include <seqan3/core/algorithm/config_base.hpp>

using namespace seqan3;

struct bar : public detail::config_base<bar>
{
public:

    friend class detail::config_access<bar>;

private:

    int state{1};
};

struct bax : public detail::config_base<bax>
{
public:
    friend class detail::config_access<bax>;

private:

    float state{2.2};
};

TEST(configurator, concept)
{
    EXPECT_TRUE((detail::configurator_concept<detail::configurator<bax, bar>>));
    // EXPECT_FALSE(configurator_concept<detail::configurator<>>);
}

TEST(configurator, tuple_size)
{
    detail::configurator<bax, bar> cfg{};

    EXPECT_EQ(2, std::tuple_size_v<decltype(cfg)>);
}

TEST(configurator, tuple_element)
{
    detail::configurator<bax, bar> cfg{};

    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, decltype(cfg)>, bax>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, decltype(cfg)>, bar>));
}

TEST(configurator, get_by_position)
{
    detail::configurator<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(get<1>(cfg), 1);
        get<1>(cfg) = 3;
        EXPECT_EQ(get<1>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<1>(cfg)), int &>));
    }

    { // const l-value
        detail::configurator<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(get<1>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<1>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configurator<bax, bar> cfg_r{cfg};
        EXPECT_EQ(get<1>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<1>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configurator<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(get<1>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<1>(std::move(cfg_rc))), int const &&>));
    }
}

TEST(configurator, get_by_position_std)
{
    detail::configurator<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(std::get<1>(cfg), 1);
        std::get<1>(cfg) = 3;
        EXPECT_EQ(std::get<1>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg)), int &>));
    }

    { // const l-value
        detail::configurator<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(std::get<1>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configurator<bax, bar> cfg_r{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configurator<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_rc))), int const &&>));
    }
}

TEST(configurator, get_by_type)
{
    detail::configurator<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(get<bar>(cfg), 1);
        get<bar>(cfg) = 3;
        EXPECT_EQ(get<bar>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(cfg)), int &>));
    }

    { // const l-value
        detail::configurator<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(get<bar>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configurator<bax, bar> cfg_r{cfg};
        EXPECT_EQ(get<bar>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configurator<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(get<bar>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(std::move(cfg_rc))), int const &&>));
    }
}

TEST(configurator, get_by_type_std)
{
    detail::configurator<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(std::get<bar>(cfg), 1);
        std::get<bar>(cfg) = 3;
        EXPECT_EQ(std::get<bar>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(cfg)), int &>));
    }

    { // const l-value
        detail::configurator<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(std::get<bar>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configurator<bax, bar> cfg_r{cfg};
        EXPECT_EQ(std::get<bar>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configurator<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(std::get<bar>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(std::move(cfg_rc))), int const &&>));
    }
}

template <size_t I>
struct foo : public detail::config_base<foo<I>>
{
public:
    friend class detail::config_access<foo<I>>;

private:

    size_t state{I};
};

TEST(configurator, replace_with)
{
    using t1 = detail::configurator<foo<0>, foo<1>, foo<2>>;
    using t2 = detail::replace_config_with_t<t1, foo<1>, foo<3>>;

    EXPECT_TRUE((std::is_same_v<t2, detail::configurator<foo<3>, foo<0>, foo<2>>>));
}

struct test_fn_impl : public detail::configurator_fn_base<test_fn_impl>
{

    template <detail::configurator_concept configurator_t>
    constexpr auto invoke(configurator_t && cfg,
                          int new_v) const
    {
        using configurator_t_ = std::remove_reference_t<configurator_t>;
        using new_cfg_t = detail::transfer_template_args_onto_t<meta::push_front<
                                                                    typename configurator_t_::type_list_type, bar>,
                                                                detail::configurator>;
        new_cfg_t new_cfg{std::forward<configurator_t>(cfg)};
        std::get<0>(new_cfg) = new_v;
        return new_cfg;
    }

    template <detail::configurator_concept configurator_t>
    constexpr auto invoke(configurator_t && cfg) const
    {
        using configurator_t_ = std::remove_reference_t<configurator_t>;
        using new_cfg_t = detail::transfer_template_args_onto_t<meta::push_front<
                                                                    typename configurator_t_::type_list_type, bar>,
                                                                detail::configurator>;
        new_cfg_t new_cfg{std::forward<configurator_t>(cfg)};
        std::get<0>(new_cfg) = 0;
        return new_cfg;
    }
};

TEST(configurator_fn, invoke_w_configurator)
{

    { // as l-value
        detail::configurator<bax> cfg;
        auto new_cfg = test_fn_impl{}(cfg, 3);

        EXPECT_EQ(get<bar>(new_cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(new_cfg), detail::configurator<bar, bax>>));
    }

    { // as r-value
        auto new_cfg = test_fn_impl{}(detail::configurator<bax>{}, 3);

        EXPECT_EQ(get<bar>(new_cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(new_cfg), detail::configurator<bar, bax>>));
    }
}

TEST(configurator_fn, pipeable_w_derived_fn)
{
    { // as l-value
        detail::configurator<bax> cfg{};
        test_fn_impl fn{};
        auto cfg_r = cfg | fn;

        EXPECT_EQ(get<0>(cfg_r), 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configurator<bar, bax>>));
    }

    { // as l-value
        test_fn_impl fn{};
        auto cfg_r = detail::configurator<bax>{} | fn;

        EXPECT_EQ(get<0>(cfg_r), 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configurator<bar, bax>>));
    }

    { // as r-value
        detail::configurator<bax> cfg{};
        auto cfg_r = cfg | test_fn_impl{};

        EXPECT_EQ(get<0>(cfg_r), 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configurator<bar, bax>>));
    }

    { // as r-value
        auto cfg_r = detail::configurator<bax>{} | test_fn_impl{};

        EXPECT_EQ(get<0>(cfg_r), 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configurator<bar, bax>>));
    }
}

TEST(configurator_fn, pipeable_w_proxy)
{
    { // as l-value
        detail::configurator<bax> cfg_{};
        int val = 3;
        auto proxy = test_fn_impl{}(val);
        auto cfg = cfg_ | proxy;

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configurator<bar, bax>>));
    }

    { // as l-value
        int val = 3;
        auto proxy = test_fn_impl{}(val);
        auto cfg = detail::configurator<bax>{} | proxy;

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configurator<bar, bax>>));
    }

    { // as r-value
        detail::configurator<bax> cfg_{};
        auto cfg = cfg_ | test_fn_impl{}(3);

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configurator<bar, bax>>));
    }

    { // as r-value
        auto cfg = detail::configurator<bax>{} | test_fn_impl{}(3);

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configurator<bar, bax>>));
    }
}
