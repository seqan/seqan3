// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// ==========================================================================

#include <gtest/gtest.h>

#include <seqan3/alphabet/all.hpp>

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

// add all alphabets here
using alphabet_types = ::testing::Types<dna4, dna5, dna15, rna4, rna5, rna15,
                                        aa27,
                                        union_composition<dna4>,
                                        union_composition<dna4, gap>,
                                        union_composition<dna5, dna5>,
                                        union_composition<dna4, dna5, gap>,
                                        /*gap,*/
                                        gapped<dna4>,
                                        gapped<dna15>,
                                        gapped<illumina18>,
                                        char, char16_t, // char32_t, too slow
                                        uint8_t, uint16_t, // uint32_t, too slow
                                        illumina18, dna4q,
                                        dot_bracket3, dssp9, wuss<>, wuss<65>,
                                        structured_rna<rna5, dot_bracket3>, structured_rna<rna4, wuss51>,
                                        structured_aa<aa27, dssp9>>;

TYPED_TEST_CASE(alphabet, alphabet_types);

TYPED_TEST(alphabet, alphabet_size)
{
    EXPECT_GT(alphabet_size_v<TypeParam>, 0u);
}

TYPED_TEST(alphabet, default_value_constructor)
{
    [[maybe_unused]] TypeParam t1;
    [[maybe_unused]] TypeParam t2{};
}

TYPED_TEST(alphabet, assign_rank)
{
    // this double checks the value initialisation
    EXPECT_EQ((assign_rank(TypeParam{}, 0)), TypeParam{});

    TypeParam t0;
    for (uint64_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
        assign_rank(t0, i);

// TODO(h-2): once we have a proper assert macro that throws instead of SIGABRTs:
//     EXPECT_THROW(assign_rank(t0, alphabet_size_v<TypeParam>));

    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(t0, 0)), TypeParam &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(TypeParam{}, 0)), TypeParam &&>));
}

TYPED_TEST(alphabet, to_rank)
{
    // this double checks the value initialisation
    EXPECT_EQ(to_rank(TypeParam{}), 0u);

    TypeParam t0;
    for (uint64_t i = 0; i < alphabet_size_v<TypeParam>; ++i)
        EXPECT_EQ((to_rank(assign_rank(t0, i))), i);

    EXPECT_TRUE((std::is_same_v<decltype(to_rank(t0)), underlying_rank_t<TypeParam>>));
}

TYPED_TEST(alphabet, copy_constructor)
{
    constexpr underlying_rank_t<TypeParam> rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    TypeParam t1;
    assign_rank(t1, rank);
    TypeParam t2{t1};
    TypeParam t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

TYPED_TEST(alphabet, move_constructor)
{
    constexpr underlying_rank_t<TypeParam> rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    TypeParam t0;
    assign_rank(t0, rank);
    TypeParam t1{t0};

    TypeParam t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    TypeParam t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

TYPED_TEST(alphabet, copy_assignment)
{
    constexpr underlying_rank_t<TypeParam> rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    TypeParam t1;
    assign_rank(t1, rank);
    TypeParam t2;
    t2 = t1;
    EXPECT_EQ(t1, t2);
}

TYPED_TEST(alphabet, move_assignment)
{
    constexpr underlying_rank_t<TypeParam> rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
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

TYPED_TEST(alphabet, swap)
{
    constexpr underlying_rank_t<TypeParam> rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    TypeParam t0;
    assign_rank(t0, rank);
    TypeParam t1{t0};
    TypeParam t2{};
    TypeParam t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

TYPED_TEST(alphabet, assign_char)
{
    using char_t = underlying_rank_t<TypeParam>;
    TypeParam t0;
    for (char_t i = std::numeric_limits<char_t>::min(); i < std::numeric_limits<char_t>::max(); ++i)
        assign_char(t0, i);

    EXPECT_TRUE((std::is_same_v<decltype(assign_char(t0, 0)), TypeParam &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char(TypeParam{}, 0)), TypeParam &&>));
}

TYPED_TEST(alphabet, to_char)
{
    TypeParam t0;
    EXPECT_TRUE((std::is_same_v<decltype(to_char(t0)), underlying_char_t<TypeParam>>));

    // more elaborate tests are done in specific alphabets
}

TYPED_TEST(alphabet, comparison_operators)
{
    if constexpr (alphabet_size_v<TypeParam> == 1)
    {
        TypeParam t0{};
        TypeParam t1{};
        EXPECT_LE(t0, t1);
        EXPECT_LE(t1, t1);
        EXPECT_EQ(t1, t1);
        EXPECT_GE(t1, t1);
        EXPECT_GE(t1, t0);
    }
    else
    {
        TypeParam t0{};
        TypeParam t1{};
        assign_rank(t0, 0);
        assign_rank(t1, 1);

        EXPECT_LT(t0, t1);
        EXPECT_LE(t0, t1);
        EXPECT_LE(t1, t1);
        EXPECT_EQ(t1, t1);
        EXPECT_NE(t0, t1);
        EXPECT_GE(t1, t1);
        EXPECT_GE(t1, t0);
        EXPECT_GT(t1, t0);
    }
}

TYPED_TEST(alphabet, concept)
{
    EXPECT_TRUE(alphabet_concept<TypeParam>);
}

#if SEQAN3_WITH_CEREAL
template <typename in_archive_t, typename out_archive_t, typename TypeParam>
void do_serialisation(TypeParam const l, std::vector<TypeParam> const & vec)
{
    // Generate unique file name.
    test::tmp_file_name filename{"alphabet_cereal_test"};
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

TYPED_TEST(alphabet, serialisation)
{
    TypeParam letter;
    if constexpr (alphabet_size_v<TypeParam> > 1)
        assign_rank(letter, 1);
    else
        assign_rank(letter, 0);

    std::vector<TypeParam> vec;
    vec.resize(10);
    if constexpr (alphabet_size_v<TypeParam> > 1)
        for (unsigned i = 0; i < 10; ++i)
            assign_rank(vec[i], i % 2);

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>(letter, vec);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(letter, vec);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>(letter, vec);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>(letter, vec);
}
#endif // SEQAN3_WITH_CEREAL

// ------------------------------------------------------------------
// constexpr tests
// ------------------------------------------------------------------

template <typename T>
class alphabet_constexpr : public ::testing::Test
{};

// add all alphabets here
using alphabet_constexpr_types = ::testing::Types<dna4, dna5, dna15, rna4, rna5, rna15,
                                                  /*aa27,*/
                                                  union_composition<dna4>,
                                                  union_composition<dna4, gap>,
                                                  union_composition<dna4, dna5, gap>,
                                                  char, char16_t, char32_t,
                                                  uint8_t, uint16_t, uint32_t,
                                                  /*gap, gapped<nucl16>, */
                                                  illumina18, dna4q,
                                                  dot_bracket3, dssp9, wuss<>, wuss<65>,
                                                  structured_rna<rna5, dot_bracket3>, structured_rna<rna4, wuss51>,
                                                  structured_aa<aa27, dssp9>>;

TYPED_TEST_CASE(alphabet_constexpr, alphabet_types);

TYPED_TEST(alphabet_constexpr, default_value_constructor)
{
    [[maybe_unused]] constexpr TypeParam t0{};
}

TYPED_TEST(alphabet_constexpr, copy_constructor)
{
    constexpr TypeParam t1{};
    constexpr TypeParam t2{t1};
    constexpr TypeParam t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

TYPED_TEST(alphabet_constexpr, move_constructor)
{
    constexpr TypeParam t0{};
    constexpr TypeParam t1{t0};

    constexpr TypeParam t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    constexpr TypeParam t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

TYPED_TEST(alphabet_constexpr, assign_rank)
{
    constexpr size_t rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    [[maybe_unused]] constexpr TypeParam t0{assign_rank(TypeParam{}, rank)};
}

TYPED_TEST(alphabet_constexpr, to_rank)
{
    constexpr size_t rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    constexpr TypeParam t0{assign_rank(TypeParam{}, rank)};
    constexpr bool b = (to_rank(t0) == rank);
    EXPECT_TRUE(b);
}

TYPED_TEST(alphabet_constexpr, copy_assignment)
{
    constexpr size_t rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    constexpr TypeParam t0{assign_rank(TypeParam{}, rank)};
    // constexpr context:
    constexpr TypeParam t3 = [&] () constexpr
    {
        TypeParam t1{assign_rank(TypeParam{}, rank)};
        TypeParam t2{};
        t2 = t1;

        return t2;
    }();
    EXPECT_EQ(t3, t0);
}

TYPED_TEST(alphabet_constexpr, move_assignment)
{
    constexpr size_t rank = (alphabet_size_v<TypeParam> == 1) ? 0 : 1;
    constexpr TypeParam t0{assign_rank(TypeParam{}, rank)};
    // constexpr context:
    constexpr TypeParam t3 = [&] () constexpr
    {
        TypeParam t1{assign_rank(TypeParam{}, rank)};
        TypeParam t2{};
        t2 = std::move(t1);

        return t2;
    }();
    EXPECT_EQ(t3, t0);
}

TYPED_TEST(alphabet_constexpr, assign_char)
{
    [[maybe_unused]] constexpr TypeParam t0{assign_char(TypeParam{}, 'A')};
}

TYPED_TEST(alphabet_constexpr, to_char)
{
    constexpr TypeParam t0{};
    [[maybe_unused]] constexpr underlying_char_t<TypeParam> c = to_char(t0);
}

TYPED_TEST(alphabet_constexpr, comparison_operators)
{
    if constexpr (alphabet_size_v<TypeParam> == 1)
    {
        constexpr TypeParam t0{};
        constexpr TypeParam t1{};
        constexpr bool b2 = (t0 <= t1);
        constexpr bool b3 = (t1 <= t1);
        constexpr bool b4 = (t1 == t1);
        constexpr bool b5 = (t1 >= t1);
        constexpr bool b6 = (t1 >= t0);

        EXPECT_TRUE(b2);
        EXPECT_TRUE(b3);
        EXPECT_TRUE(b4);
        EXPECT_TRUE(b5);
        EXPECT_TRUE(b6);
    }
    else
    {
        constexpr TypeParam t0{assign_rank(TypeParam{}, 0)};
        constexpr TypeParam t1{assign_rank(TypeParam{}, 1)};
        constexpr bool b1 = (t0 <  t1);
        constexpr bool b2 = (t0 <= t1);
        constexpr bool b3 = (t1 <= t1);
        constexpr bool b4 = (t1 == t1);
        constexpr bool b5 = (t1 >= t1);
        constexpr bool b6 = (t1 >= t0);
        constexpr bool b7 = (t1 >  t0);
        constexpr bool b8 = (t0 != t1);

        EXPECT_TRUE(b1);
        EXPECT_TRUE(b2);
        EXPECT_TRUE(b3);
        EXPECT_TRUE(b4);
        EXPECT_TRUE(b5);
        EXPECT_TRUE(b6);
        EXPECT_TRUE(b7);
        EXPECT_TRUE(b8);
    }
}
