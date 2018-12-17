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

#include <functional>
#include <random>
#include <type_traits>

#include <seqan3/core/metafunction/function.hpp>
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

using namespace seqan3;
using namespace seqan3::align_cfg;

template <typename value_t>
struct dummy_gap : seq_end_gap_base<value_t>
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

using static_end_gap_types = ::testing::Types<first_seq_leading<std::true_type>,
                                              first_seq_leading<std::false_type>,
                                              first_seq_trailing<std::true_type>,
                                              first_seq_trailing<std::false_type>,
                                              second_seq_leading<std::true_type>,
                                              second_seq_leading<std::false_type>,
                                              second_seq_trailing<std::true_type>,
                                              second_seq_trailing<std::false_type>>;

TYPED_TEST_CASE(static_end_gap_test, static_end_gap_types);

template <typename type>
class dynamic_end_gap_test : public ::testing::Test
{};

using dynamic_end_gap_types = ::testing::Types<first_seq_leading<bool>,
                                               first_seq_trailing<bool>,
                                               second_seq_leading<bool>,
                                               second_seq_trailing<bool>>;

TYPED_TEST_CASE(dynamic_end_gap_test, dynamic_end_gap_types);

TEST(seq_end_gap_base, aggregate)
{
    EXPECT_TRUE((std::is_aggregate_v<dummy_gap<std::true_type>>));
    EXPECT_TRUE((std::is_aggregate_v<dummy_gap<bool>>));
}

TEST(first_seq_leading, deduction)
{
    { // static
        first_seq_leading tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), first_seq_leading<std::true_type>>), true);
        first_seq_leading tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), first_seq_leading<std::false_type>>), true);
    }

    { // dynamic
        first_seq_leading tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), first_seq_leading<bool>>), true);
        first_seq_leading tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), first_seq_leading<bool>>), true);
    }
}

TEST(first_seq_trailing, deduction)
{
    { // static
        first_seq_trailing tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), first_seq_trailing<std::true_type>>), true);
        first_seq_trailing tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), first_seq_trailing<std::false_type>>), true);
    }

    { // dynamic
        first_seq_trailing tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), first_seq_trailing<bool>>), true);
        first_seq_trailing tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), first_seq_trailing<bool>>), true);
    }
}

TEST(second_seq_leading, deduction)
{
    { // static
        second_seq_leading tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), second_seq_leading<std::true_type>>), true);
        second_seq_leading tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), second_seq_leading<std::false_type>>), true);
    }

    { // dynamic
        second_seq_leading tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), second_seq_leading<bool>>), true);
        second_seq_leading tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), second_seq_leading<bool>>), true);
    }
}

TEST(second_seq_trailing, deduction)
{
    { // static
        second_seq_trailing tmp{std::true_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp), second_seq_trailing<std::true_type>>), true);
        second_seq_trailing tmp2{std::false_type{}};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), second_seq_trailing<std::false_type>>), true);
    }

    { // dynamic
        second_seq_trailing tmp{true};
        EXPECT_EQ((std::is_same_v<decltype(tmp), second_seq_trailing<bool>>), true);
        second_seq_trailing tmp2{false};
        EXPECT_EQ((std::is_same_v<decltype(tmp2), second_seq_trailing<bool>>), true);
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
        EXPECT_TRUE((std::is_nothrow_default_constructible_v<end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_copy_constructible_v<end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_move_constructible_v<end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_copy_assignable_v<end_gaps<>>));
        EXPECT_TRUE((std::is_nothrow_move_assignable_v<end_gaps<>>));
    }

    { // one element
        EXPECT_TRUE((std::is_nothrow_default_constructible_v<end_gaps<first_seq_leading<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_copy_constructible_v<end_gaps<first_seq_leading<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_move_constructible_v<end_gaps<first_seq_leading<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_copy_assignable_v<end_gaps<first_seq_leading<std::true_type>>>));
        EXPECT_TRUE((std::is_nothrow_move_assignable_v<end_gaps<first_seq_leading<std::true_type>>>));

        EXPECT_TRUE((std::is_nothrow_constructible_v<end_gaps<first_seq_leading<std::true_type>>>));
    }

    { // four elements
        EXPECT_TRUE((std::is_nothrow_default_constructible_v<end_gaps<first_seq_leading<std::true_type>,
                                                                      second_seq_leading<bool>,
                                                                      first_seq_trailing<std::false_type>,
                                                                      second_seq_trailing<bool>>>));
        EXPECT_TRUE((std::is_nothrow_copy_constructible_v<end_gaps<first_seq_leading<std::true_type>,
                                                                   second_seq_leading<bool>,
                                                                   first_seq_trailing<std::false_type>,
                                                                   second_seq_trailing<bool>>>));
        EXPECT_TRUE((std::is_nothrow_move_constructible_v<end_gaps<first_seq_leading<std::true_type>,
                                                                   second_seq_leading<bool>,
                                                                   first_seq_trailing<std::false_type>,
                                                                   second_seq_trailing<bool>>>));
        EXPECT_TRUE((std::is_nothrow_copy_assignable_v<end_gaps<first_seq_leading<std::true_type>,
                                                                second_seq_leading<bool>,
                                                                first_seq_trailing<std::false_type>,
                                                                second_seq_trailing<bool>>>));
        EXPECT_TRUE((std::is_nothrow_move_assignable_v<end_gaps<first_seq_leading<std::true_type>,
                                                                second_seq_leading<bool>,
                                                                first_seq_trailing<std::false_type>,
                                                                second_seq_trailing<bool>>>));
        EXPECT_TRUE((std::is_nothrow_constructible_v<end_gaps<first_seq_leading<std::true_type>,
                                                              second_seq_leading<bool>,
                                                              first_seq_trailing<std::false_type>,
                                                              second_seq_trailing<bool>>>));
    }
}

TEST(end_gaps, deduction)
{
    { // default
        end_gaps eg{};
        EXPECT_EQ((std::is_same_v<decltype(eg), end_gaps<>>), true);
    }

    { // one element
        end_gaps eg{second_seq_trailing{true}};
        using foo = end_gaps<second_seq_trailing<bool>>;
        EXPECT_EQ((std::is_same_v<decltype(eg), foo>), true);
    }

    { // multiple elements
        end_gaps eg{second_seq_trailing{true}, first_seq_leading{std::true_type{}}, second_seq_leading{std::false_type{}}};
        using foo = end_gaps<second_seq_trailing<bool>, first_seq_leading<std::true_type>, second_seq_leading<std::false_type>>;
        EXPECT_EQ((std::is_same_v<decltype(eg), foo>), true);
    }
}

TEST(end_gaps, access)
{
    { // default
        end_gaps eg{};

        EXPECT_EQ(eg[0], false);
        EXPECT_EQ(eg[1], false);
        EXPECT_EQ(eg[2], false);
        EXPECT_EQ(eg[3], false);
    }

    { // custom
        end_gaps eg{second_seq_trailing{true}, first_seq_leading{std::true_type{}}, second_seq_leading{std::false_type{}}};

        EXPECT_EQ(eg[0], true);
        EXPECT_EQ(eg[1], false);
        EXPECT_EQ(eg[2], false);
        EXPECT_EQ(eg[3], true);
    }
}

TEST(end_gaps, static_query)
{
    { // default
        end_gaps eg{};

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
        end_gaps eg{second_seq_trailing{true}, first_seq_leading{std::true_type{}}, second_seq_leading{std::false_type{}}};

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
    end_gaps eg{second_seq_trailing{true}, first_seq_leading{std::true_type{}}, second_seq_leading{std::false_type{}}};

    constexpr bool seq1_l = decltype(eg)::get_static<0>();
    constexpr bool seq2_l = decltype(eg)::get_static<2>();

    EXPECT_EQ(seq1_l, true);
    EXPECT_EQ(eg[0], true);
    EXPECT_EQ(seq2_l, false);
    EXPECT_EQ(eg[2], false);
}

TEST(end_gaps, all_ends_free)
{
    using test = end_gaps<first_seq_leading<std::true_type>,
                          first_seq_trailing<std::true_type>,
                          second_seq_leading<std::true_type>,
                          second_seq_trailing<std::true_type>>;
    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(all_ends_free)>, test>), true);
}

TEST(end_gaps, none_ends_free)
{
    using test = end_gaps<first_seq_leading<std::false_type>,
                          first_seq_trailing<std::false_type>,
                          second_seq_leading<std::false_type>,
                          second_seq_trailing<std::false_type>>;

    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(none_ends_free)>, test>), true);
}

TEST(end_gaps, seq1_ends_free)
{
    using test = end_gaps<first_seq_leading<std::true_type>,
                          first_seq_trailing<std::true_type>,
                          second_seq_leading<std::false_type>,
                          second_seq_trailing<std::false_type>>;
    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(seq1_ends_free)>, test>), true);
}

TEST(end_gaps, seq2_ends_free)
{
    using test = end_gaps<first_seq_leading<std::false_type>,
                          first_seq_trailing<std::false_type>,
                          second_seq_leading<std::true_type>,
                          second_seq_trailing<std::true_type>>;
    EXPECT_EQ((std::is_same_v<std::remove_const_t<decltype(seq2_ends_free)>, test>), true);
}

TEST(align_cfg_aligned_ends, is_aggregate)
{
    EXPECT_TRUE((std::is_aggregate_v<aligned_ends<end_gaps<>>>));
}

TEST(align_cfg_aligned_ends, id)
{
    align_cfg::aligned_ends cfg{all_ends_free};

    EXPECT_EQ(static_cast<uint8_t>(decltype(cfg)::id),
              static_cast<uint8_t>(detail::align_config_id::aligned_ends));
}

TEST(align_cfg_aligned_ends, value)
{
    align_cfg::aligned_ends cfg{seq1_ends_free};
    using type = decltype(cfg.value);

    using test = end_gaps<first_seq_leading<std::true_type>,
                          first_seq_trailing<std::true_type>,
                          second_seq_leading<std::false_type>,
                          second_seq_trailing<std::false_type>>;

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
        align_cfg::aligned_ends elem{all_ends_free};
        configuration cfg{elem};

        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[0]), true);
        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[1]), true);
        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[2]), true);
        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[3]), true);
    }

    {
        configuration cfg{align_cfg::aligned_ends{all_ends_free}};

        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[0]), true);
        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[1]), true);
        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[2]), true);
        EXPECT_EQ((get<align_cfg::aligned_ends>(cfg).value[3]), true);
    }
}
