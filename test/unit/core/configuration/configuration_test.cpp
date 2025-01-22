// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>

#include "configuration_mock.hpp"

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
    seqan3::configuration cfg0{};      // empty
    seqan3::configuration cfg1{bax{}}; // one element

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
        EXPECT_TRUE((std::is_same_v<decltype(std::get<1>(std::move(cfg_rc))), bar const &&>));
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
        EXPECT_TRUE((std::is_same_v<decltype(std::get<bar>(std::move(cfg_rc))), bar const &&>));
    }
}

TEST(configuration, get_by_type_template)
{
    seqan3::configuration cfg = bar{1} | foobar<>{std::vector{0, 1, 2, 3}};

    { // l-value
        EXPECT_RANGE_EQ(seqan3::get<foobar>(cfg).value, (std::vector{0, 1, 2, 3}));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(cfg)), foobar<> &>));
    }

    { // const l-value
        seqan3::configuration<bar, foobar<>> const cfg_c{cfg};
        EXPECT_RANGE_EQ(seqan3::get<foobar>(cfg_c).value, (std::vector{0, 1, 2, 3}));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(cfg_c)), foobar<> const &>));
    }

    { // r-value
        seqan3::configuration<bar, foobar<>> cfg_r{cfg};
        EXPECT_RANGE_EQ(seqan3::get<foobar>(std::move(cfg_r)).value, (std::vector{0, 1, 2, 3}));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(std::move(cfg_r))), foobar<> &&>));
    }

    { // const r-value
        seqan3::configuration<bar, foobar<>> const cfg_cr{cfg};
        EXPECT_RANGE_EQ(seqan3::get<foobar>(std::move(cfg_cr)).value, (std::vector{0, 1, 2, 3}));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::get<foobar>(std::move(cfg_cr))), foobar<> const &&>));
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

TEST(configuration, append_configuration_element)
{
    {
        seqan3::configuration<foo, bar> cfg{};
        auto new_cfg = cfg.append(bax{});
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<foo, bar, bax>>));
    }

    {
        seqan3::configuration<foo, bar> cfg{};
        bax b{};
        auto new_cfg = cfg.append(b);
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<foo, bar, bax>>));
    }

    {
        seqan3::configuration<foo, bar> cfg{};
        bax b{};
        auto new_cfg = cfg.append(std::as_const(b));
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<foo, bar, bax>>));
    }

    {
        seqan3::configuration<> cfg{};
        bax b{};
        auto new_cfg = cfg.append(std::as_const(b));
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<bax>>));
    }
}

TEST(configuration, append_configuration)
{
    {
        seqan3::configuration<foo, bar> cfg{};
        auto new_cfg = cfg.append(seqan3::configuration{bax{}});
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<foo, bar, bax>>));
    }

    {
        seqan3::configuration<foo, bar> cfg{};
        seqan3::configuration cfg2{bax{}};
        auto new_cfg = cfg.append(cfg2);
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<foo, bar, bax>>));
    }

    {
        seqan3::configuration<foo, bar> cfg{};
        seqan3::configuration cfg2{bax{}};
        auto new_cfg = cfg.append(std::as_const(cfg2));
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<foo, bar, bax>>));
    }

    {
        seqan3::configuration<> cfg{};
        seqan3::configuration cfg2{bax{}};
        auto new_cfg = cfg.append(std::as_const(cfg2));
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<bax>>));
    }

    {
        seqan3::configuration<> cfg{};
        auto new_cfg = cfg.append(seqan3::configuration<>{});
        EXPECT_TRUE((std::same_as<decltype(new_cfg), seqan3::configuration<>>));
    }
}

TEST(configuration, remove_by_type)
{
    {
        seqan3::configuration<foo, bax, bar> cfg{};

        // remove middle
        EXPECT_SAME_TYPE(decltype(cfg.template remove<bax>()), (seqan3::configuration<foo, bar>));
        // remove head
        EXPECT_SAME_TYPE(decltype(cfg.template remove<bar>()), (seqan3::configuration<foo, bax>));
        // remove tail
        EXPECT_SAME_TYPE(decltype(cfg.template remove<foo>()), (seqan3::configuration<bax, bar>));

        seqan3::configuration<foo> single_cfg{};
        EXPECT_SAME_TYPE(decltype(single_cfg.template remove<foo>()), seqan3::configuration<>);
    }

    { // const cfg
        seqan3::configuration<foo, bax, bar> const cfg{};

        // remove middle
        EXPECT_SAME_TYPE(decltype(cfg.template remove<bax>()), (seqan3::configuration<foo, bar>));
        // remove head
        EXPECT_SAME_TYPE(decltype(cfg.template remove<bar>()), (seqan3::configuration<foo, bax>));
        // remove tail
        EXPECT_SAME_TYPE(decltype(cfg.template remove<foo>()), (seqan3::configuration<bax, bar>));

        seqan3::configuration<foo> const single_cfg{};
        EXPECT_SAME_TYPE(decltype(single_cfg.template remove<foo>()), seqan3::configuration<>);
    }
}

TEST(configuration, remove_by_type_template)
{
    {
        seqan3::configuration<foo, foobar<>, bar> cfg{};

        EXPECT_SAME_TYPE(decltype(cfg.template remove<foobar>()), (seqan3::configuration<foo, bar>));

        seqan3::configuration<foobar<>> single_cfg{};
        EXPECT_SAME_TYPE(decltype(single_cfg.template remove<foobar>()), seqan3::configuration<>);
    }

    { //const cfg
        seqan3::configuration<foo, foobar<>, bar> const cfg{};

        EXPECT_SAME_TYPE(decltype(cfg.template remove<foobar>()), (seqan3::configuration<foo, bar>));

        seqan3::configuration<foobar<>> const single_cfg{};
        EXPECT_SAME_TYPE(decltype(single_cfg.template remove<foobar>()), seqan3::configuration<>);
    }
}

TEST(configuration, get_or_by_type)
{
    seqan3::configuration cfg = bax{2.2} | bar{1};

    { // l-value
        EXPECT_FLOAT_EQ(cfg.get_or(bax{1.3}).value, 2.2);
        EXPECT_SAME_TYPE(decltype(cfg.get_or(bax{1.3})), bax &);
        EXPECT_EQ(cfg.get_or(foo{"test"}).value, "test");
    }

    { // const l-value
        EXPECT_FLOAT_EQ(std::as_const(cfg).get_or(bax{1.3}).value, 2.2);
        EXPECT_SAME_TYPE(decltype(std::as_const(cfg).get_or(bax{1.3})), bax const &);
        EXPECT_EQ(std::as_const(cfg).get_or(foo{"test"}).value, "test");
    }

    { // r-value
        EXPECT_FLOAT_EQ(std::move(seqan3::configuration<bax, bar>{cfg}).get_or(bax{1.3}).value, 2.2);
        EXPECT_SAME_TYPE(decltype(seqan3::configuration<bax, bar>{cfg}.get_or(bax{1.3})), bax &&);
        EXPECT_EQ(std::move(seqan3::configuration<bax, bar>{cfg}).get_or(foo{"test"}).value, "test");
    }

    { // const r-value
        seqan3::configuration<bax, bar> const cfg_cr{cfg};
        EXPECT_FLOAT_EQ(std::move(cfg_cr).get_or(bax{1.3}).value, 2.2);
        seqan3::configuration<bax, bar> const cfg_cr2{cfg};
        EXPECT_SAME_TYPE(decltype(std::move(cfg_cr2).get_or(bax{1.3})), bax const &&);
        seqan3::configuration<bax, bar> const cfg_cr3{cfg};
        EXPECT_EQ(std::move(cfg_cr3).get_or(foo{"test"}).value, "test");
    }
}

TEST(configuration, get_or_by_type_constexpr)
{
    constexpr seqan3::configuration cfg = bax{2.2} | bar{1};

    { // l-value
        constexpr auto element = cfg.get_or(bax{1.3});
        EXPECT_FLOAT_EQ(element.value, 2.2);
    }

    { // const l-value
        constexpr auto element = std::as_const(cfg).get_or(bax{1.3});
        EXPECT_FLOAT_EQ(element.value, 2.2);
    }

    { // r-value
        constexpr auto element = seqan3::configuration{cfg}.get_or(bax{1.3});
        EXPECT_FLOAT_EQ(element.value, 2.2);
    }

    { // const r-value
        constexpr seqan3::configuration const cfg_c{cfg};
        constexpr auto element = std::move(cfg_c).get_or(bax{1.3});
        EXPECT_FLOAT_EQ(element.value, 2.2);
    }
}

TEST(configuration, get_or_by_type_template)
{
    seqan3::configuration cfg = bar{1} | foobar<>{std::vector{0, 1, 2, 3}};

    using double_vec_t = std::vector<double>;
    using alternative_t = foobar<double_vec_t>;
    alternative_t alternative{double_vec_t{3.3}};

    { // l-value
        EXPECT_RANGE_EQ(cfg.get_or(alternative_t{double_vec_t{3.3}}).value, (std::vector{0, 1, 2, 3}));
        EXPECT_RANGE_EQ(cfg.get_or(alternative).value, (std::vector{0, 1, 2, 3}));
        EXPECT_EQ(cfg.get_or(foo{"test"}).value, "test");
        EXPECT_SAME_TYPE(decltype(cfg.get_or(alternative)), foobar<std::vector<int>> &);
    }

    { // const l-value
        EXPECT_RANGE_EQ(std::as_const(cfg).get_or(alternative_t{double_vec_t{3.3}}).value, (std::vector{0, 1, 2, 3}));
        EXPECT_RANGE_EQ(std::as_const(cfg).get_or(alternative).value, (std::vector{0, 1, 2, 3}));
        EXPECT_EQ(std::as_const(cfg).get_or(foo{"test"}).value, "test");
        EXPECT_SAME_TYPE(decltype(std::as_const(cfg).get_or(alternative)), foobar<std::vector<int>> const &);
    }

    { // r-value;
        EXPECT_RANGE_EQ((seqan3::configuration<bar, foobar<>>{cfg}.get_or(alternative_t{double_vec_t{3.3}}).value),
                        (std::vector{0, 1, 2, 3}));
        EXPECT_RANGE_EQ((seqan3::configuration<bar, foobar<>>{cfg}.get_or(alternative).value),
                        (std::vector{0, 1, 2, 3}));
        EXPECT_EQ((seqan3::configuration<bar, foobar<>>{cfg}.get_or(foo{"test"}).value), "test");
        EXPECT_TRUE((std::same_as<decltype(seqan3::configuration<bar, foobar<>>{cfg}.get_or(alternative)),
                                  foobar<std::vector<int>> &&>));
    }

    { // const r-value
        seqan3::configuration<bar, foobar<>> const cfg_cr{cfg};
        EXPECT_RANGE_EQ(std::move(cfg_cr).get_or(alternative_t{double_vec_t{3.3}}).value, (std::vector{0, 1, 2, 3}));
        seqan3::configuration<bar, foobar<>> const cfg_cr2{cfg};
        EXPECT_RANGE_EQ(std::move(cfg_cr2).get_or(alternative).value, (std::vector{0, 1, 2, 3}));
        seqan3::configuration<bar, foobar<>> const cfg_cr3{cfg};
        EXPECT_EQ(std::move(cfg_cr3).get_or(foo{"test"}).value, "test");
        seqan3::configuration<bar, foobar<>> const cfg_cr4{cfg};
        EXPECT_TRUE(
            (std::same_as<decltype(std::move(cfg_cr4).get_or(alternative)), foobar<std::vector<int>> const &&>));
    }
}

TEST(configuration, get_or_by_template_type_constexpr)
{
    constexpr seqan3::configuration cfg = bar{1} | foobar<bool>{true};

    { // l-value
        constexpr auto element = cfg.get_or(foobar<int>{1});
        EXPECT_EQ(element.value, true);
    }

    { // const l-value
        constexpr auto element = std::as_const(cfg).get_or(foobar<int>{1});
        EXPECT_EQ(element.value, true);
    }

    { // r-value
        constexpr auto element = seqan3::configuration{cfg}.get_or(foobar<int>{1});
        EXPECT_EQ(element.value, true);
    }

    { // const r-value
        constexpr seqan3::configuration const cfg_c{cfg};
        constexpr auto element = std::move(cfg_c).get_or(foobar<int>{1});
        EXPECT_EQ(element.value, true);
    }
}

TEST(configuration, get_or_perfectly_forwarded_alternative)
{
    seqan3::configuration<bar, foo> cfg{};

    using alternative_t = foobar<std::vector<double>>;
    alternative_t alternative{std::vector<double>{3.3}};
    alternative_t const const_alternative{alternative};

    EXPECT_RANGE_EQ(cfg.get_or(alternative_t{std::vector<double>{3.3}}).value, (std::vector<double>{3.3}));
    EXPECT_SAME_TYPE(decltype(cfg.get_or(alternative)), alternative_t &);
    EXPECT_SAME_TYPE(decltype(cfg.get_or(const_alternative)), alternative_t const &);
    EXPECT_SAME_TYPE(decltype(cfg.get_or(alternative_t{std::vector<double>{3.3}})), alternative_t);
    EXPECT_SAME_TYPE(decltype(cfg.get_or(std::move(const_alternative))), alternative_t const);
}
