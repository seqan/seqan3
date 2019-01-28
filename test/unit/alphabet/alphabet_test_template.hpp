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
#include <seqan3/test/pretty_printing.hpp>

#if SEQAN3_WITH_CEREAL
#include <seqan3/test/tmp_filename.hpp>

#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#endif // SEQAN3_WITH_CEREAL

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
    EXPECT_EQ((assign_rank(TypeParam{}, 0)), TypeParam{});

    TypeParam t0;
    for (uint64_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
        assign_rank(t0, i);

// TODO(h-2): once we have a proper assert macro that throws instead of SIGABRTs:
//     EXPECT_THROW(assign_rank(t0, alphabet_size_v<TypeParam>));

    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(t0, 0)), TypeParam &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(TypeParam{}, 0)), TypeParam>));
}

TYPED_TEST_P(alphabet, global_to_rank)
{
    // this double checks the value initialisation
    EXPECT_EQ(to_rank(TypeParam{}), 0u);

    TypeParam t0;
    for (uint64_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
        EXPECT_EQ((to_rank(assign_rank(t0, i))), i);

    EXPECT_TRUE((std::is_same_v<decltype(to_rank(t0)), underlying_rank_t<TypeParam>>));
}

TYPED_TEST_P(alphabet, copy_constructor)
{
    // the module operation ensures that the result is within the valid rank range;
    // it will be in the most cases 1 except for alphabets like seqan3::gap where it will be 0
    constexpr underlying_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t1;
    assign_rank(t1, rank);
    TypeParam t2{t1};
    TypeParam t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

TYPED_TEST_P(alphabet, move_constructor)
{
    constexpr underlying_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t0;
    assign_rank(t0, rank);
    TypeParam t1{t0};

    TypeParam t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    TypeParam t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(alphabet, copy_assignment)
{
    constexpr underlying_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t1;
    assign_rank(t1, rank);
    TypeParam t2;
    t2 = t1;
    EXPECT_EQ(t1, t2);
}

TYPED_TEST_P(alphabet, move_assignment)
{
    constexpr underlying_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t0;
    assign_rank(t0, rank);
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
    constexpr underlying_rank_t<TypeParam> rank = 1 % alphabet_size_v<TypeParam>;
    TypeParam t0;
    assign_rank(t0, rank);
    TypeParam t1{t0};
    TypeParam t2{};
    TypeParam t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

TYPED_TEST_P(alphabet, global_assign_char)
{
    using char_t = underlying_char_t<TypeParam>;
    char_t i = std::numeric_limits<char_t>::min();
    char_t j = std::numeric_limits<char_t>::max();

    TypeParam t0;
    for (; i < j; ++i)
        assign_char(t0, i);

    EXPECT_TRUE((std::is_same_v<decltype(assign_char(t0, 0)), TypeParam &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char(TypeParam{}, 0)), TypeParam>));
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
    for (underlying_char_t<TypeParam> c :
         ranges::view::iota(ptrdiff_t{std::numeric_limits<underlying_char_t<TypeParam>>::min()},
                            ptrdiff_t{std::numeric_limits<underlying_char_t<TypeParam>>::max()} + 1))
    {
        if (char_is_valid_for<TypeParam>(c))
            EXPECT_NO_THROW(assign_char_strict(TypeParam{}, c));
        else
            EXPECT_THROW(assign_char_strict(TypeParam{}, c), invalid_char_assignment);
    }
}

TYPED_TEST_P(alphabet, global_to_char)
{
    TypeParam t0;
    EXPECT_TRUE((std::is_same_v<decltype(to_char(t0)), underlying_char_t<TypeParam>>));

    // more elaborate tests are done in specific alphabets

}

TYPED_TEST_P(alphabet, comparison_operators)
{
    TypeParam t0{};
    TypeParam t1{};

    assign_rank(t0, 0);
    assign_rank(t1, 1 % alphabet_size_v<TypeParam>);

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
    EXPECT_TRUE(alphabet_concept<TypeParam>);
    EXPECT_TRUE(alphabet_concept<TypeParam &>);
    EXPECT_TRUE(alphabet_concept<TypeParam &&>);
}

TYPED_TEST_P(alphabet, debug_streaming)
{
    std::ostringstream o;
    debug_stream.set_underlying_stream(o);

    debug_stream << TypeParam{};

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
            for (uint64_t i = 0; i < alphabet_size_v<TypeParam>/2; ++i)
            {
                assign_rank(t0, i);
                ASSERT_EQ(h(t0), i);
            }
        }
        else
        {
            for (uint64_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
            {
                assign_rank(t0, i);
                ASSERT_EQ(h(t0), i);
            }
        }
    }
    {
        std::vector<TypeParam> text;
        text.reserve(4);
        for (uint64_t i = 0; i < 4; ++i)
        {
            text.push_back(assign_rank(TypeParam{}, 0));
        }
        std::hash<decltype(text)> h{};
        ASSERT_EQ(h(text), 0u);
    }
}

#if SEQAN3_WITH_CEREAL
template <typename in_archive_t, typename out_archive_t, typename TypeParam>
void do_serialisation(TypeParam const l, std::vector<TypeParam> const & vec)
{
    // Generate unique file name.
    test::tmp_filename filename{"alphabet_cereal_test"};
    {
        std::ofstream os{filename.get_path(), std::ios::binary};
        out_archive_t oarchive{os};
        oarchive(l);
        oarchive(vec);
    }

    {
        TypeParam in_l{};
        std::vector<TypeParam> in_vec;

        std::ifstream is{filename.get_path(), std::ios::binary};
        in_archive_t iarchive{is};
        iarchive(in_l);
        iarchive(in_vec);
        EXPECT_EQ(l, in_l);
        EXPECT_EQ(vec, in_vec);
    }
}
#endif // SEQAN3_WITH_CEREAL

TYPED_TEST_P(alphabet, serialisation)
{
#if SEQAN3_WITH_CEREAL
    TypeParam letter;

    assign_rank(letter, 1 % alphabet_size_v<TypeParam>);

    std::vector<TypeParam> vec;
    vec.resize(10);
    for (unsigned i = 0; i < 10; ++i)
        assign_rank(vec[i], i % alphabet_size_v<TypeParam>);

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>(letter, vec);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(letter, vec);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>(letter, vec);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>(letter, vec);
#endif // SEQAN3_WITH_CEREAL
}

REGISTER_TYPED_TEST_CASE_P(alphabet, alphabet_size, default_value_constructor, global_assign_rank, global_to_rank,
    copy_constructor, move_constructor, copy_assignment, move_assignment, swap, global_assign_char, global_to_char,
    global_char_is_valid_for, global_assign_char_strict, comparison_operators, concept_check, debug_streaming, hash,
    serialisation);
