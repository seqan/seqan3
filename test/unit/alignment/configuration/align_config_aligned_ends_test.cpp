// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <functional>
#include <random>
#include <type_traits>

#include <seqan3/core/type_traits/function.hpp>
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

template <typename value_t>
struct dummy_gap : seqan3::sequence_end_gap_specifier_base<value_t>
{};

template <typename type>
class static_end_gap_test : public ::testing::Test
{
public:

    static constexpr bool expected_value = static_end_gap_test::determine_value(type{});

private:

    template <template <typename> typename gap_element,
              typename value_t>
    static constexpr bool determine_value(gap_element<value_t>)
    {
        return value_t::value;
    }
};

using static_end_gap_types = ::testing::Types<seqan3::front_end_first<std::true_type>,
                                              seqan3::front_end_first<std::false_type>,
                                              seqan3::back_end_first<std::true_type>,
                                              seqan3::back_end_first<std::false_type>,
                                              seqan3::front_end_second<std::true_type>,
                                              seqan3::front_end_second<std::false_type>,
                                              seqan3::back_end_second<std::true_type>,
                                              seqan3::back_end_second<std::false_type>>;

TYPED_TEST_SUITE(static_end_gap_test, static_end_gap_types, );

template <typename type>
class dynamic_end_gap_test : public ::testing::Test
{};

using dynamic_end_gap_types = ::testing::Types<seqan3::front_end_first<bool>,
                                               seqan3::back_end_first<bool>,
                                               seqan3::front_end_second<bool>,
                                               seqan3::back_end_second<bool>>;

TYPED_TEST_SUITE(dynamic_end_gap_test, dynamic_end_gap_types, );

TEST(sequence_end_gap_specifier_base, aggregate)
{
    EXPECT_TRUE((std::is_aggregate_v<dummy_gap<std::true_type>>));
    EXPECT_TRUE((std::is_aggregate_v<dummy_gap<bool>>));
}

TEST(front_end_first, deduction)
{
    { // static
        seqan3::front_end_first tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::front_end_first<std::true_type>>), true);
        seqan3::front_end_first tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::front_end_first<std::false_type>>), true);
    }

    { // dynamic
        seqan3::front_end_first tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::front_end_first<bool>>), true);
        seqan3::front_end_first tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::front_end_first<bool>>), true);
    }
}

TEST(back_end_first, deduction)
{
    { // static
        seqan3::back_end_first tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::back_end_first<std::true_type>>), true);
        seqan3::back_end_first tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::back_end_first<std::false_type>>), true);
    }

    { // dynamic
        seqan3::back_end_first tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::back_end_first<bool>>), true);
        seqan3::back_end_first tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::back_end_first<bool>>), true);
    }
}

TEST(front_end_second, deduction)
{
    { // static
        seqan3::front_end_second tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::front_end_second<std::true_type>>), true);
        seqan3::front_end_second tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::front_end_second<std::false_type>>), true);
    }

    { // dynamic
        seqan3::front_end_second tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::front_end_second<bool>>), true);
        seqan3::front_end_second tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::front_end_second<bool>>), true);
    }
}

TEST(back_end_second, deduction)
{
    { // static
        seqan3::back_end_second tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::back_end_second<std::true_type>>), true);
        seqan3::back_end_second tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::back_end_second<std::false_type>>), true);
    }

    { // dynamic
        seqan3::back_end_second tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), seqan3::back_end_second<bool>>), true);
        seqan3::back_end_second tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), seqan3::back_end_second<bool>>), true);
    }
}

TYPED_TEST(static_end_gap_test, aggreagte)
{
    EXPECT_TRUE((std::is_aggregate_v<TypeParam>));
}

TYPED_TEST(static_end_gap_test, get_value)
{
    TypeParam instance{};
    EXPECT_EQ(instance(), TestFixture::expected_value);
}

TYPED_TEST(dynamic_end_gap_test, construction)
{
    EXPECT_TRUE((std::is_aggregate_v<TypeParam>));
}

TYPED_TEST(dynamic_end_gap_test, get_value)
{
    { // true
        TypeParam instance{true};
        EXPECT_EQ(instance(), true);
    }
    { // false
        TypeParam instance{false};
        EXPECT_EQ(instance(), false);
    }
}

TEST(end_gaps, construction)
{
    { // empty
        EXPECT_TRUE((std::is_nothrow_default_constructible_v<seqan3::end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_copy_constructible_v<seqan3::end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_move_constructible_v<seqan3::end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_copy_assignable_v<seqan3::end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_move_assignable_v<seqan3::end_gaps<>>));
    }

    { // one element
        EXPECT_TRUE((std::is_nothrow_default_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_copy_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_move_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_copy_assignable_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_move_assignable_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>>>));

        EXPECT_TRUE((std::is_nothrow_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>>>));
    }

    { // four elements
        EXPECT_TRUE((std::is_nothrow_default_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                                                              seqan3::front_end_second<bool>,
                                                                              seqan3::back_end_first<std::false_type>,
                                                                              seqan3::back_end_second<bool>>>));
        EXPECT_TRUE((std::is_nothrow_copy_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                                                           seqan3::front_end_second<bool>,
                                                                           seqan3::back_end_first<std::false_type>,
                                                                           seqan3::back_end_second<bool>>>));
        EXPECT_TRUE((std::is_nothrow_move_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                                                           seqan3::front_end_second<bool>,
                                                                           seqan3::back_end_first<std::false_type>,
                                                                           seqan3::back_end_second<bool>>>));
        EXPECT_TRUE((std::is_nothrow_copy_assignable_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                                                        seqan3::front_end_second<bool>,
                                                                        seqan3::back_end_first<std::false_type>,
                                                                        seqan3::back_end_second<bool>>>));
        EXPECT_TRUE((std::is_nothrow_move_assignable_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                                                        seqan3::front_end_second<bool>,
                                                                        seqan3::back_end_first<std::false_type>,
                                                                        seqan3::back_end_second<bool>>>));
        EXPECT_TRUE((std::is_nothrow_constructible_v<seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                                                      seqan3::front_end_second<bool>,
                                                                      seqan3::back_end_first<std::false_type>,
                                                                      seqan3::back_end_second<bool>>>));
    }

    { // from lvalue
        seqan3::front_end_first fsl{std::true_type{}};
        EXPECT_TRUE(seqan3::end_gaps{fsl}[0]);
    }
}

TEST(end_gaps, deduction)
{
    { // default
        seqan3::end_gaps eg{};
        EXPECT_EQ((std::is_same_v<decltype(eg), seqan3::end_gaps<>>), true);
    }

    { // one element
        seqan3::end_gaps eg{seqan3::back_end_second{true}};
        using foo = seqan3::end_gaps<seqan3::back_end_second<bool>>;
        EXPECT_EQ((std::is_same_v<decltype(eg), foo>), true);
    }

    { // multiple elements
        seqan3::end_gaps eg{seqan3::back_end_second{true},
                            seqan3::front_end_first{std::true_type{}},
                            seqan3::front_end_second{std::false_type{}}};
        using foo = seqan3::end_gaps<seqan3::back_end_second<bool>,
                                     seqan3::front_end_first<std::true_type>,
                                     seqan3::front_end_second<std::false_type>>;
        EXPECT_EQ((std::is_same_v<decltype(eg), foo>), true);
    }
}

TEST(end_gaps, access)
{
    { // default
        seqan3::end_gaps eg{};

        EXPECT_EQ(eg[0], false);
        EXPECT_EQ(eg[1], false);
        EXPECT_EQ(eg[2], false);
        EXPECT_EQ(eg[3], false);
    }

    { // custom
        seqan3::end_gaps eg{seqan3::back_end_second{true},
                            seqan3::front_end_first{std::true_type{}},
                            seqan3::front_end_second{std::false_type{}}};

        EXPECT_EQ(eg[0], true);
        EXPECT_EQ(eg[1], false);
        EXPECT_EQ(eg[2], false);
        EXPECT_EQ(eg[3], true);
    }
}

TEST(end_gaps, static_query)
{
    { // default
        seqan3::end_gaps eg{};

        constexpr bool seq1_l = decltype(eg)::is_static<0>();
        constexpr bool seq1_t = decltype(eg)::is_static<1>();
        constexpr bool seq2_l = decltype(eg)::is_static<2>();
        constexpr bool seq2_t = decltype(eg)::is_static<3>();
        EXPECT_EQ(seq1_l, false);
        EXPECT_EQ(seq1_t, false);
        EXPECT_EQ(seq2_l, false);
        EXPECT_EQ(seq2_t, false);
    }

    { // custom
        seqan3::end_gaps eg{seqan3::back_end_second{true},
                            seqan3::front_end_first{std::true_type{}},
                            seqan3::front_end_second{std::false_type{}}};

        constexpr bool seq1_l = decltype(eg)::is_static<0>();
        constexpr bool seq1_t = decltype(eg)::is_static<1>();
        constexpr bool seq2_l = decltype(eg)::is_static<2>();
        constexpr bool seq2_t = decltype(eg)::is_static<3>();

        EXPECT_EQ(seq1_l, true);
        EXPECT_EQ(seq1_t, false);
        EXPECT_EQ(seq2_l, true);
        EXPECT_EQ(seq2_t, false);
    }
}

TEST(end_gaps, static_access)
{
    seqan3::end_gaps eg{seqan3::back_end_second{true},
                        seqan3::front_end_first{std::true_type{}},
                        seqan3::front_end_second{std::false_type{}}};

    constexpr bool seq1_l = decltype(eg)::get_static<0>();
    constexpr bool seq2_l = decltype(eg)::get_static<2>();

    EXPECT_EQ(seq1_l, true);
    EXPECT_EQ(eg[0], true);
    EXPECT_EQ(seq2_l, false);
    EXPECT_EQ(eg[2], false);
}

TEST(end_gaps, free_ends_all)
{
    using test = seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                  seqan3::back_end_first<std::true_type>,
                                  seqan3::front_end_second<std::true_type>,
                                  seqan3::back_end_second<std::true_type>>;
    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(seqan3::free_ends_all)>, test>), true);
}

TEST(end_gaps, free_ends_none)
{
    using test = seqan3::end_gaps<seqan3::front_end_first<std::false_type>,
                                  seqan3::back_end_first<std::false_type>,
                                  seqan3::front_end_second<std::false_type>,
                                  seqan3::back_end_second<std::false_type>>;

    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(seqan3::free_ends_none)>, test>), true);
}

TEST(end_gaps, free_ends_first)
{
    using test = seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                  seqan3::back_end_first<std::true_type>,
                                  seqan3::front_end_second<std::false_type>,
                                  seqan3::back_end_second<std::false_type>>;
    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(seqan3::free_ends_first)>, test>), true);
}

TEST(end_gaps, free_ends_second)
{
    using test = seqan3::end_gaps<seqan3::front_end_first<std::false_type>,
                                  seqan3::back_end_first<std::false_type>,
                                  seqan3::front_end_second<std::true_type>,
                                  seqan3::back_end_second<std::true_type>>;
    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(seqan3::free_ends_second)>, test>), true);
}

TEST(align_cfg_aligned_ends, is_aggregate)
{
    EXPECT_TRUE((std::is_aggregate_v<seqan3::align_cfg::aligned_ends<seqan3::end_gaps<>>>));
}

TEST(align_cfg_aligned_ends, id)
{
    seqan3::align_cfg::aligned_ends cfg{seqan3::free_ends_all};

    EXPECT_EQ(static_cast<uint8_t>(decltype(cfg)::id),
              static_cast<uint8_t>(seqan3::detail::align_config_id::aligned_ends));
}

TEST(align_cfg_aligned_ends, value)
{
    seqan3::align_cfg::aligned_ends cfg{seqan3::free_ends_first};
    using type = decltype(cfg.value);

    using test = seqan3::end_gaps<seqan3::front_end_first<std::true_type>,
                                  seqan3::back_end_first<std::true_type>,
                                  seqan3::front_end_second<std::false_type>,
                                  seqan3::back_end_second<std::false_type>>;

    EXPECT_EQ((std::is_same_v<type, test>), true);

    EXPECT_EQ(cfg.value[0], true);
    EXPECT_EQ(cfg.value[1], true);
    EXPECT_EQ(cfg.value[2], false);
    EXPECT_EQ(cfg.value[3], false);

    EXPECT_EQ(type::is_static<0>(), true);
    EXPECT_EQ(type::get_static<0>(), true);
    EXPECT_EQ(type::is_static<2>(), true);
    EXPECT_EQ(type::get_static<2>(), false);
}

TEST(align_cfg_aligned_ends, configuration)
{
    {
        seqan3::align_cfg::aligned_ends elem{seqan3::free_ends_all};
        seqan3::configuration cfg{elem};

        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[0]), true);
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[1]), true);
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[2]), true);
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[3]), true);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::aligned_ends{seqan3::free_ends_all}};

        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[0]), true);
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[1]), true);
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[2]), true);
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::aligned_ends>(cfg).value[3]), true);
    }
}
