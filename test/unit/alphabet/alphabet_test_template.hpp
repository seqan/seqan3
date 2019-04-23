// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/iota.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/exception.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/test/cereal.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

template <typename T>
class alphabet : public ::testing::Test
{};

TYPED_TEST_CASE_P(alphabet);

TYPED_TEST_P(alphabet, alphabet_size)
{
    EXPECT_GT(alphabet_size_v<TypeParam>, 0u);
}

TYPED_TEST_P(alphabet, default_value_constructor)
{
    [[maybe_unused]] TypeParam t1;
    [[maybe_unused]] TypeParam t2{};
}

TYPED_TEST_P(alphabet, global_assign_rank)
{
    // this double checks the value initialisation
    EXPECT_EQ((assign_rank_to(0, TypeParam{})), TypeParam{});

    TypeParam t0;
    for (uint64_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
        assign_rank_to(i, t0);

// TODO(h-2): once we have a proper assert macro that throws instead of SIGABRTs:
//     EXPECT_THROW(assign_rank_to(alphabet_size_v<TypeParam>, t0));

    EXPECT_TRUE((std::is_same_v<decltype(assign_rank_to(0, t0)), TypeParam &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank_to(0, TypeParam{})), TypeParam>));
}

TYPED_TEST_P(alphabet, global_to_rank)
{
    // this double checks the value initialisation
    EXPECT_EQ(to_rank(TypeParam{}), 0u);

    TypeParam t0;
    for (uint64_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
        EXPECT_EQ((to_rank(assign_rank_to(i, t0))), i);

    EXPECT_TRUE((std::is_same_v<decltype(to_rank(t0)), alphabet_rank_t<TypeParam>>));
}

TYPED_TEST_P(alphabet, copy_constructor)
{
    // the module operation ensures that the result is within the valid rank range;
    // it will be in the most cases 1 except for alphabets like seqan3::gap where it will be 0
    constexpr alphabet_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t1;
    assign_rank_to(rank, t1);
    TypeParam t2{t1};
    TypeParam t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

TYPED_TEST_P(alphabet, move_constructor)
{
    constexpr alphabet_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t0;
    assign_rank_to(rank, t0);
    TypeParam t1{t0};

    TypeParam t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    TypeParam t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(alphabet, copy_assignment)
{
    constexpr alphabet_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t1;
    assign_rank_to(rank, t1);
    TypeParam t2;
    t2 = t1;
    EXPECT_EQ(t1, t2);
}

TYPED_TEST_P(alphabet, move_assignment)
{
    constexpr alphabet_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t0;
    assign_rank_to(rank, t0);
    TypeParam t1{t0};
    TypeParam t2;
    TypeParam t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(alphabet, swap)
{
    constexpr alphabet_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t0;
    assign_rank_to(rank, t0);
    TypeParam t1{t0};
    TypeParam t2{};
    TypeParam t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

TYPED_TEST_P(alphabet, global_assign_char)
{
    using char_t = alphabet_char_t<TypeParam>;
    char_t i = std::numeric_limits<char_t>::min();
    char_t j = std::numeric_limits<char_t>::max();

    TypeParam t0;
    for (; i < j; ++i)
        assign_char_to(i, t0);

    EXPECT_TRUE((std::is_same_v<decltype(assign_char_to(0, t0)), TypeParam &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_to(0, TypeParam{})), TypeParam>));
}

TYPED_TEST_P(alphabet, global_char_is_valid_for) // only test negative example for most; more inside specialised tests
{
    if constexpr (alphabet_size_v<TypeParam> < 255) // includes most of our alphabets, but not the adaptations!
    {
        EXPECT_FALSE((char_is_valid_for<TypeParam>(0))); // for none of our alphabets char{0} is valid
    }
}

TYPED_TEST_P(alphabet, global_assign_char_strict)
{
    for (alphabet_char_t<TypeParam> c :
         std::view::iota(ptrdiff_t{std::numeric_limits<alphabet_char_t<TypeParam>>::min()},
                            ptrdiff_t{std::numeric_limits<alphabet_char_t<TypeParam>>::max()} + 1))
    {
        if (char_is_valid_for<TypeParam>(c))
            EXPECT_NO_THROW(assign_char_strictly_to(c, TypeParam{}));
        else
            EXPECT_THROW(assign_char_strictly_to(c, TypeParam{}), invalid_char_assignment);
    }
}

TYPED_TEST_P(alphabet, global_to_char)
{
    TypeParam t0;
    EXPECT_TRUE((std::is_same_v<decltype(to_char(t0)), alphabet_char_t<TypeParam>>));

    // more elaborate tests are done in specific alphabets

}

TYPED_TEST_P(alphabet, comparison_operators)
{
    TypeParam t0{};
    TypeParam t1{};

    assign_rank_to(0, t0);
    assign_rank_to(1 % alphabet_size_v<TypeParam>, t1);

    EXPECT_EQ(t0, t0);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t1, t0);

    if constexpr (alphabet_size_v<TypeParam> == 1)
    {
        EXPECT_EQ(t0, t1);
    }
    else
    {
        EXPECT_LT(t0, t1);
        EXPECT_NE(t0, t1);
        EXPECT_GT(t1, t0);
    }
}

TYPED_TEST_P(alphabet, concept_check)
{
    EXPECT_TRUE(Alphabet<TypeParam>);
    EXPECT_TRUE(Alphabet<TypeParam &>);
    EXPECT_TRUE(Alphabet<TypeParam &&>);
}

TYPED_TEST_P(alphabet, debug_streaming)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    my_stream << TypeParam{};

    o.flush();
    EXPECT_EQ(o.str().size(), 1u);
}

TYPED_TEST_P(alphabet, hash)
{
    {
        TypeParam t0{};
        std::hash<TypeParam> h{};
        if constexpr (std::Same<TypeParam, char>)
        {
            for (size_t i = 0; i < alphabet_size_v<TypeParam>/2; ++i)
            {
                assign_rank_to(i, t0);
                ASSERT_EQ(h(t0), i);
            }
        }
        else
        {
            for (size_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
            {
                assign_rank_to(i, t0);
                ASSERT_EQ(h(t0), i);
            }
        }
    }
    {
        std::vector<TypeParam> text;
        text.reserve(4);
        for (size_t i = 0; i < 4; ++i)
        {
            text.push_back(assign_rank_to(0, TypeParam{}));
        }
        std::hash<decltype(text)> h{};
        ASSERT_EQ(h(text), 0u);
    }
    {
        std::hash<TypeParam const> h{};
        if constexpr (std::Same<TypeParam, char>)
        {
            for (size_t i = 0; i < alphabet_size_v<TypeParam>/2; ++i)
            {
                TypeParam const t0 = assign_rank_to(i, TypeParam{});
                ASSERT_EQ(h(t0), i);
            }
        }
        else
        {
            for (size_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
            {
                TypeParam const t0 = assign_rank_to(i, TypeParam{});
                ASSERT_EQ(h(t0), i);
            }
        }
    }
    {
        std::vector<TypeParam> const text(4, assign_rank_to(0, TypeParam{}));
        std::hash<decltype(text)> h{};
        ASSERT_EQ(h(text), 0u);
    }
}

TYPED_TEST_P(alphabet, serialisation)
{
    TypeParam letter;

    assign_rank_to(1 % alphabet_size_v<TypeParam>, letter);
    test::do_serialisation(letter);

    std::vector<TypeParam> vec;
    vec.resize(10);
    for (unsigned i = 0; i < 10; ++i)
        assign_rank_to(i % alphabet_size_v<TypeParam>, vec[i]);
    test::do_serialisation(vec);
}

REGISTER_TYPED_TEST_CASE_P(alphabet, alphabet_size, default_value_constructor, global_assign_rank, global_to_rank,
    copy_constructor, move_constructor, copy_assignment, move_assignment, swap, global_assign_char, global_to_char,
    global_char_is_valid_for, global_assign_char_strict, comparison_operators, concept_check, debug_streaming, hash,
    serialisation);
