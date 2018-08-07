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
#include <seqan3/core/algorithm/config_element_base.hpp>

using namespace seqan3;

struct bar : public detail::config_element_base<bar>
{
public:

    friend class detail::config_element_access<bar>;

    bar() : detail::config_element_base<bar>{}
    {
        ++default_counter;
    }

    bar(bar const & b) : detail::config_element_base<bar>{b}, state(b.state)
    {
        ++copy_counter;
    }

    bar(bar && b) : detail::config_element_base<bar>{std::move(b)}, state(b.state)
    {
        ++move_counter;
    }

    bar & operator=(bar const & b)
    {
        state = b.state;
        ++copy_counter;
        return *this;
    }

    bar & operator=(bar && b)
    {
        state = std::move(b.state);
        ++move_counter;
        return *this;
    }

    ~bar() = default;

    static inline size_t default_counter{0};
    static inline size_t copy_counter{0};
    static inline size_t move_counter{0};

    static void reset_counter()
    {
        default_counter = 0;
        copy_counter = 0;
        move_counter = 0;
    }

private:

    int state{1};
};

struct bax : public detail::config_element_base<bax>
{
public:
    friend class detail::config_element_access<bax>;

    bax() : detail::config_element_base<bax>{}
    {
        ++default_counter;
    }

    bax(bax const & b) : detail::config_element_base<bax>{b}, state(b.state)
    {
        ++copy_counter;
    }

    bax(bax && b) : detail::config_element_base<bax>{std::move(b)}, state(b.state)
    {
        ++move_counter;
    }

    bax & operator=(bax const & b)
    {
        state = b.state;
        ++copy_counter;
        return *this;
    }

    bax & operator=(bax && b)
    {
        state = std::move(b.state);
        ++move_counter;
        return *this;
    }

    ~bax() = default;

    static inline size_t default_counter{0};
    static inline size_t copy_counter{0};
    static inline size_t move_counter{0};

    static void reset_counter()
    {
        default_counter = 0;
        copy_counter = 0;
        move_counter = 0;
    }

private:

    float state{2.2};
};

TEST(configuration, metafunction)
{
    EXPECT_TRUE((detail::is_algorithm_configuration<detail::configuration<bax, bar>>::value));
    EXPECT_TRUE((detail::is_algorithm_configuration_v<detail::configuration<bax, bar>>));
    EXPECT_FALSE((detail::is_algorithm_configuration<type_list<bax>>::value));
    EXPECT_FALSE((detail::is_algorithm_configuration_v<type_list<bax>>));
}

TEST(configuration, tuple_size)
{
    detail::configuration<bax, bar> cfg{};

    EXPECT_EQ(2, std::tuple_size_v<decltype(cfg)>);
}

TEST(configuration, tuple_element)
{
    detail::configuration<bax, bar> cfg{};

    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, decltype(cfg)>, bax>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, decltype(cfg)>, bar>));
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
        EXPECT_EQ(get<1>(cfg), 1);
        get<1>(cfg) = 3;
        EXPECT_EQ(get<1>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<1>(cfg)), int &>));
    }

    { // const l-value
        detail::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(get<1>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<1>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(get<1>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<1>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(get<1>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<1>(std::move(cfg_rc))), int const &&>));
    }
}

TEST(configuration, get_by_position_std)
{
    detail::configuration<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(std::get<1>(cfg), 1);
        std::get<1>(cfg) = 3;
        EXPECT_EQ(std::get<1>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg)), int &>));
    }

    { // const l-value
        detail::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(std::get<1>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(std::get<1>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_rc))), int const &&>));
    }
}

TEST(configuration, get_by_type)
{
    detail::configuration<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(get<bar>(cfg), 1);
        get<bar>(cfg) = 3;
        EXPECT_EQ(get<bar>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(cfg)), int &>));
    }

    { // const l-value
        detail::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(get<bar>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(get<bar>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(get<bar>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(get<bar>(std::move(cfg_rc))), int const &&>));
    }
}

TEST(configuration, get_by_type_std)
{
    detail::configuration<bax, bar> cfg{};

    { // l-value
        EXPECT_EQ(std::get<bar>(cfg), 1);
        std::get<bar>(cfg) = 3;
        EXPECT_EQ(std::get<bar>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(cfg)), int &>));
    }

    { // const l-value
        detail::configuration<bax, bar> const cfg_c{cfg};
        EXPECT_EQ(std::get<bar>(cfg_c), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(cfg_c)), int const &>));
    }

    { // r-value
        detail::configuration<bax, bar> cfg_r{cfg};
        EXPECT_EQ(std::get<bar>(std::move(cfg_r)), 3);
        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(std::move(cfg_r))), int &&>));
    }

    { // const r-value
        detail::configuration<bax, bar> const cfg_rc{cfg};
        EXPECT_EQ(std::get<bar>(std::move(cfg_rc)), 3);

        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(std::move(cfg_rc))), int const &&>));
    }
}

TEST(configuration, construction_w_lvalue_cfg)
{
    bar::reset_counter();
    bax::reset_counter();

    detail::configuration<bar> cfg;

    EXPECT_EQ(bar::default_counter, 1);
    EXPECT_EQ(bar::copy_counter, 0);
    EXPECT_EQ(bar::move_counter, 0);
    EXPECT_EQ(bax::default_counter, 0);
    EXPECT_EQ(bax::copy_counter, 0);
    EXPECT_EQ(bax::move_counter, 0);

    detail::configuration<bax, bar> cfg2{cfg};

    EXPECT_EQ(bar::default_counter, 2);
    EXPECT_EQ(bar::copy_counter, 1);
    EXPECT_EQ(bar::move_counter, 0);
    EXPECT_EQ(bax::default_counter, 1);
    EXPECT_EQ(bax::copy_counter, 0);
    EXPECT_EQ(bax::move_counter, 0);

    EXPECT_EQ(get<1>(cfg2), 1);
}

TEST(configuration, construction_w_rvalue_cfg)
{
    bar::reset_counter();
    bax::reset_counter();

    detail::configuration<bar> cfg;

    EXPECT_EQ(bar::default_counter, 1);
    EXPECT_EQ(bar::copy_counter, 0);
    EXPECT_EQ(bar::move_counter, 0);
    EXPECT_EQ(bax::default_counter, 0);
    EXPECT_EQ(bax::copy_counter, 0);
    EXPECT_EQ(bax::move_counter, 0);

    detail::configuration<bax, bar> cfg2{std::move(cfg)};

    EXPECT_EQ(bar::default_counter, 2);
    EXPECT_EQ(bar::copy_counter, 0);
    EXPECT_EQ(bar::move_counter, 1);
    EXPECT_EQ(bax::default_counter, 1);
    EXPECT_EQ(bax::copy_counter, 0);
    EXPECT_EQ(bax::move_counter, 0);

    EXPECT_EQ(get<1>(cfg2), 1);
}

template <size_t I>
struct foo : public detail::config_element_base<foo<I>>
{
public:
    friend class detail::config_element_access<foo<I>>;

private:

    size_t state{I};
};

TEST(configuration, replace_with)
{
    using t1 = detail::configuration<foo<0>, foo<1>, foo<2>>;
    using t2 = detail::replace_config_with_t<t1, foo<1>, foo<3>>;

    EXPECT_TRUE((std::is_same_v<t2, detail::configuration<foo<3>, foo<0>, foo<2>>>));
}

struct bar_fn_impl : public detail::configuration_fn_base<bar_fn_impl>
{

    // using detail::configuration_fn_base<bar_fn_impl>::configuration_fn_base;

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg,
                          int new_v) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        using configuration_t_ = std::remove_reference_t<configuration_t>;
        using new_cfg_t = detail::push_front_config_t<configuration_t_, bar>;

        new_cfg_t new_cfg{std::forward<configuration_t>(cfg)};
        std::get<0>(new_cfg) = new_v;
        return new_cfg;
    }

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        using configuration_t_ = std::remove_reference_t<configuration_t>;
        using new_cfg_t = detail::push_front_config_t<configuration_t_, bar>;

        new_cfg_t new_cfg{std::forward<configuration_t>(cfg)};
        std::get<0>(new_cfg) = 0;
        return new_cfg;
    }
};

struct bax_fn_impl : public detail::configuration_fn_base<bax_fn_impl>
{
    // using detail::configuration_fn_base<bax_fn_impl>::configuration_fn_base;

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg,
                          float new_v) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        using configuration_t_ = std::remove_reference_t<configuration_t>;
        using new_cfg_t = detail::push_front_config_t<configuration_t_, bax>;

        new_cfg_t new_cfg{std::forward<configuration_t>(cfg)};
        std::get<0>(new_cfg) = new_v;
        return new_cfg;
    }

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        using configuration_t_ = std::remove_reference_t<configuration_t>;
        using new_cfg_t = detail::push_front_config_t<configuration_t_, bax>;

        new_cfg_t new_cfg{std::forward<configuration_t>(cfg)};
        std::get<0>(new_cfg) = 0.0;
        return new_cfg;
    }
};

TEST(configuration, template_deduction_from_proxy)
{
    detail::configuration cfg = bar_fn_impl{}(3);

    EXPECT_EQ(get<0>(cfg), 3);
    EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar>>));
}

TEST(configuration, template_deduction_from_variable)
{
    detail::configuration cfg = bar_fn_impl{};

    EXPECT_EQ(get<0>(cfg), 0);
    EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar>>));
}

TEST(configuration_fn, invoke_w_configuration)
{

    { // as l-value
        detail::configuration<bax> cfg;
        auto new_cfg = bar_fn_impl{}(cfg, 3);

        EXPECT_EQ(get<bar>(new_cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(new_cfg), detail::configuration<bar, bax>>));
    }

    { // as r-value
        auto new_cfg = bar_fn_impl{}(detail::configuration<bax>{}, 3);

        EXPECT_EQ(get<bar>(new_cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(new_cfg), detail::configuration<bar, bax>>));
    }
}

TEST(configuration_fn, pipeable_w_derived_fn)
{
    { // as l-value
        detail::configuration<bax> cfg{};
        bar_fn_impl fn{};
        auto cfg_r = cfg | fn;

        EXPECT_EQ(get<0>(cfg_r), 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configuration<bar, bax>>));
    }

    { // as l-value
        bar_fn_impl fn{};
        auto cfg_r = detail::configuration<bax>{} | fn;

        EXPECT_EQ(get<0>(cfg_r), 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configuration<bar, bax>>));
    }

    { // as r-value
        detail::configuration<bax> cfg{};
        auto cfg_r = cfg | bar_fn_impl{};

        EXPECT_EQ(get<0>(cfg_r), 0);
        EXPECT_TRUE((std::is_same_v<decltype(cfg_r), detail::configuration<bar, bax>>));
    }

    { // as r-value
        auto cfg_r = detail::configuration<bax>{} | bar_fn_impl{};

        EXPECT_EQ(get<0>(cfg_r), 0);
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

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }

    { // as l-value
        int val = 3;
        auto proxy = bar_fn_impl{}(val);
        auto cfg = detail::configuration<bax>{} | proxy;

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }

    { // as r-value
        detail::configuration<bax> cfg_{};
        auto cfg = cfg_ | bar_fn_impl{}(3);

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }

    { // as r-value
        auto cfg = detail::configuration<bax>{} | bar_fn_impl{}(3);

        EXPECT_EQ(get<0>(cfg), 3);
        EXPECT_TRUE((std::is_same_v<decltype(cfg), detail::configuration<bar, bax>>));
    }
}

TEST(configuration_fn, pipable_fn_fn)
{
    auto cfg = bar_fn_impl{} | bax_fn_impl{};

    EXPECT_EQ(get<bar>(cfg), 0);
    EXPECT_FLOAT_EQ(get<bax>(cfg), 0.0);
}

TEST(configuration_fn, pipable_fn_proxy)
{
    auto cfg = bar_fn_impl{} | bax_fn_impl{}(3.0);

    EXPECT_EQ(get<bar>(cfg), 0);
    EXPECT_FLOAT_EQ(get<bax>(cfg), 3.0);
}

TEST(configuration_fn, pipable_proxy_fn)
{
    auto cfg = bar_fn_impl{}(2) | bax_fn_impl{};

    EXPECT_EQ(get<bar>(cfg), 2);
    EXPECT_FLOAT_EQ(get<bax>(cfg), 0.0);
}

TEST(configuration_fn, pipable_proxy_proxy)
{
    auto cfg = bar_fn_impl{}(2) | bax_fn_impl{}(3.0);

    EXPECT_EQ(get<bar>(cfg), 2);
    EXPECT_FLOAT_EQ(get<bax>(cfg), 3.0);
}
