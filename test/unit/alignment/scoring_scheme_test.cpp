// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/core/type_list.hpp>

#if SEQAN3_WITH_CEREAL
#include <seqan3/test/tmp_filename.hpp>

#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#endif // SEQAN3_WITH_CEREAL

using namespace seqan3;

template <typename T>
class generic : public ::testing::Test
{
public:
    using alph_t = std::conditional_t<detail::is_type_specialisation_of_v<T, nucleotide_scoring_scheme>,
                                      dna15,
                                      aa27>;
};

using scoring_scheme_types = ::testing::Types<nucleotide_scoring_scheme<>,
                                              nucleotide_scoring_scheme<float>,
                                              aminoacid_scoring_scheme<>,
                                              aminoacid_scoring_scheme<int>>;

TYPED_TEST_CASE(generic, scoring_scheme_types);


TEST(nucleotide_scoring_scheme, template_argument_deduction)
{

    {
        nucleotide_scoring_scheme scheme;
        EXPECT_TRUE((std::is_same_v<decltype(scheme), nucleotide_scoring_scheme<int8_t>>));
    }

    {
        nucleotide_scoring_scheme scheme{};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), nucleotide_scoring_scheme<int8_t>>));
    }

    {
        nucleotide_scoring_scheme scheme{match_score{6}, mismatch_score{-4}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), nucleotide_scoring_scheme<int8_t>>));
    }

    {
        std::array<std::array<int16_t, 15>, 15> m;
        nucleotide_scoring_scheme scheme{m};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), nucleotide_scoring_scheme<int16_t>>));
    }
}

TEST(aminoacid_scoring_scheme, template_argument_deduction)
{
    {
        aminoacid_scoring_scheme scheme;
        EXPECT_TRUE((std::is_same_v<decltype(scheme), aminoacid_scoring_scheme<int8_t>>));
    }

    {
        aminoacid_scoring_scheme scheme{};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), aminoacid_scoring_scheme<int8_t>>));
    }


    {
        aminoacid_scoring_scheme scheme{match_score{6}, mismatch_score{-4}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), aminoacid_scoring_scheme<int8_t>>));
    }

    {
        std::array<std::array<int16_t, 27>, 27> m;
        aminoacid_scoring_scheme scheme{m};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), aminoacid_scoring_scheme<int16_t>>));
    }

    {
        aminoacid_scoring_scheme scheme{aminoacid_similarity_matrix::BLOSUM62};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), aminoacid_scoring_scheme<int8_t>>));
    }
}

// ------------------------------------------------------------------
// generic test
// ------------------------------------------------------------------

TYPED_TEST(generic, concept_check)
{
    using alph_t = typename TestFixture::alph_t;
    EXPECT_TRUE((scoring_scheme_concept<TypeParam, alph_t>));
    EXPECT_TRUE((scoring_scheme_concept<TypeParam const, alph_t>));
    EXPECT_TRUE((scoring_scheme_concept<TypeParam const &, alph_t>));
}

TYPED_TEST(generic, member_types)
{
    using alph_t = typename TestFixture::alph_t;

    TypeParam scheme{};

    using score_t = typename TypeParam::score_type;
    using matrix_size_t = typename TypeParam::matrix_size_type;
    using matrix_t = typename TypeParam::matrix_type;
    constexpr auto matrix_size = TypeParam::matrix_size;

    if constexpr (std::is_same_v<TypeParam, nucleotide_scoring_scheme<float>>)
    {
        EXPECT_TRUE((std::is_same_v<score_t, float>));
    }
    else if constexpr (std::is_same_v<TypeParam, aminoacid_scoring_scheme<int>>)
    {
        EXPECT_TRUE((std::is_same_v<score_t, int>));
    }
    else
    {
        EXPECT_TRUE((std::is_same_v<score_t, int8_t>));
    }

    EXPECT_EQ(matrix_size, alphabet_size_v<alph_t>);
    EXPECT_TRUE((std::is_same_v<std::remove_const_t<decltype(matrix_size)>, matrix_size_t>));
    EXPECT_TRUE((std::is_same_v<matrix_size_t, uint8_t>));

    EXPECT_TRUE((std::is_same_v<matrix_t, std::array<std::array<score_t, matrix_size>, matrix_size>>));
}

TYPED_TEST(generic, simple_score)
{
    using alph_t = typename TestFixture::alph_t;

    // Test constructor
    TypeParam scheme{match_score{6}, mismatch_score{-4}};
    // Test set function
    scheme.set_simple_scheme(match_score{5}, mismatch_score{-3});

    for (uint8_t i = 0; i < alphabet_size_v<alph_t>; ++i)
    {
        for (uint8_t j = 0; j < alphabet_size_v<alph_t>; ++j)
        {
            int8_t expected = i == j ? 5 : -3;
            EXPECT_EQ(expected, scheme.score(assign_rank(alph_t{}, i), assign_rank(alph_t{}, j)));
        }
    }
}

TYPED_TEST(generic, simple_score_failure)
{
    if constexpr (std::is_same_v<typename TypeParam::score_type, int8_t>)
    {
        // Test constructor
        EXPECT_THROW((TypeParam{match_score{600}, mismatch_score{-4}}),
                     std::invalid_argument);

        TypeParam scheme{};
        // Test set function
        EXPECT_THROW((scheme.set_simple_scheme(match_score{-150}, mismatch_score{-3})),
                     std::invalid_argument);
    }
    else
    {
        // Test constructor
        EXPECT_NO_THROW((TypeParam{match_score{600}, mismatch_score{-4}}));

        TypeParam scheme{};
        // Test set function
        EXPECT_NO_THROW((scheme.set_simple_scheme(match_score{-150}, mismatch_score{-3})));
    }
}

TYPED_TEST(generic, hamming)
{
    using alph_t = typename TestFixture::alph_t;

    // Test constructor
    TypeParam scheme{};
    // Test set function
    scheme.set_hamming_distance();

    for (uint8_t i = 0; i < alphabet_size_v<alph_t>; ++i)
    {
        for (uint8_t j = 0; j < alphabet_size_v<alph_t>; ++j)
        {
            int8_t expected = i == j ? 0 : -1;
            EXPECT_EQ(expected, scheme.score(assign_rank(alph_t{}, i), assign_rank(alph_t{}, j)));
        }
    }
}

TYPED_TEST(generic, custom)
{
    using alph_t = typename TestFixture::alph_t;

    typename TypeParam::matrix_type matrix;

    for (uint8_t i = 0; i < alphabet_size_v<alph_t>; ++i)
        for (uint8_t j = 0; j < alphabet_size_v<alph_t>; ++j)
            matrix[i][j] = i * i + j;

    // Test constructor
    TypeParam scheme{matrix};
    // Test set function
    scheme.set_custom_matrix(matrix);

    if constexpr (detail::is_type_specialisation_of_v<TypeParam, aminoacid_scoring_scheme>)
    {
        EXPECT_EQ(0*0+0,    scheme.score(aa27::A, aa27::A));
        EXPECT_EQ(0*0+2,    scheme.score(aa27::A, aa27::C));
        EXPECT_EQ(2*2+0,    scheme.score(aa27::C, aa27::A));
        EXPECT_EQ(8*8+8,    scheme.score(aa27::I, aa27::I));
        EXPECT_EQ(0*0+13,   scheme.score(aa27::A, aa27::N));
        EXPECT_EQ(2*2+1,    scheme.score(aa27::C, aa27::B));
    } else
    {
        EXPECT_EQ(0*0+0,    scheme.score(dna15::A, dna15::A)); // A is 0th
        EXPECT_EQ(0*0+2,    scheme.score(dna15::A, dna15::C)); // C is 2nd
        EXPECT_EQ(2*2+0,    scheme.score(dna15::C, dna15::A));
        EXPECT_EQ(3*3+3,    scheme.score(dna15::D, dna15::D)); // D is 3rd
        EXPECT_EQ(8*8+0,    scheme.score(dna15::N, dna15::A)); // N is 8th
        EXPECT_EQ(0*0+8,    scheme.score(dna15::A, dna15::N));
    }
}

TYPED_TEST(generic, convertability)
{
    using alph_t = typename TestFixture::alph_t;

    typename TypeParam::matrix_type matrix;

    for (uint8_t i = 0; i < alphabet_size_v<alph_t>; ++i)
        for (uint8_t j = 0; j < alphabet_size_v<alph_t>; ++j)
            matrix[i][j] = i * i + j;

    TypeParam scheme{};
    scheme.set_custom_matrix(matrix);

    if constexpr (detail::is_type_specialisation_of_v<TypeParam, aminoacid_scoring_scheme>)
    {
        using aa_types = type_list<aa27, aa20>;
        meta::for_each(aa_types{}, [&] (auto && aa) constexpr
        {
            using nucl_t = std::decay_t<decltype(aa)>;

            EXPECT_EQ(scheme.score(aa27::C,  aa27::G),
                      scheme.score(nucl_t::C, nucl_t::G));
            EXPECT_EQ(scheme.score(aa27::T,  nucl_t::T),
                      scheme.score(nucl_t::T, aa27::T));
        });
    } else
    {

        using nucl_types = type_list<dna4, dna5, dna15, rna4, rna5, rna15>;
        meta::for_each(nucl_types{}, [&] (auto && nucl) constexpr
        {
            using nucl_t = std::decay_t<decltype(nucl)>;

            EXPECT_EQ(scheme.score(dna15::C,  dna15::G),
                      scheme.score(nucl_t::C, nucl_t::G));
            EXPECT_EQ(scheme.score(dna15::T,  dna15::T),
                      scheme.score(nucl_t::T, nucl_t::T));
            EXPECT_EQ(scheme.score(dna15::A,  dna15::C),
                      scheme.score(nucl_t::A, nucl_t::C));

            EXPECT_EQ(scheme.score(dna15::C,  nucl_t::G),
                      scheme.score(nucl_t::C, dna15::G));
            EXPECT_EQ(scheme.score(dna15::C,  nucl_t::A),
                      scheme.score(nucl_t::C, dna15::A));
        });

    }
}

#if SEQAN3_WITH_CEREAL
template <typename in_archive_t, typename out_archive_t, typename TypeParam>
void do_serialisation(TypeParam const l)
{
    // This makes sure the file is also deleted if an exception is thrown in one of the tests below
    // Generate unique file name.
    test::tmp_filename filename{"scoring_scheme_cereal_test"};

    {
        std::ofstream os{filename.get_path().string(), std::ios::binary};
        out_archive_t oarchive{os};
        oarchive(l);
    }

    {
        TypeParam in_l{};
        std::ifstream is{filename.get_path().string(), std::ios::binary};
        in_archive_t iarchive{is};
        iarchive(in_l);
        EXPECT_EQ(l, in_l);
    }
}

TYPED_TEST(generic, serialisation)
{
    TypeParam scheme1;
    scheme1.set_hamming_distance();

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);

    scheme1.set_simple_scheme(match_score{11}, mismatch_score{-7});

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);
}
#endif // SEQAN3_WITH_CEREAL

// ------------------------------------------------------------------
// aminoacid test
// ------------------------------------------------------------------

template <typename T>
class aminoacid : public ::testing::Test {};

using aa_scheme_types = ::testing::Types<aminoacid_scoring_scheme<>,
                                         aminoacid_scoring_scheme<int>>;

TYPED_TEST_CASE(aminoacid, aa_scheme_types);

TYPED_TEST(aminoacid, similarity_matrix)
{
    // Test constructor
    aminoacid_scoring_scheme scheme{aminoacid_similarity_matrix::BLOSUM30};
    EXPECT_EQ( 4,    scheme.score(aa27::A, aa27::A));
    EXPECT_EQ(-3,    scheme.score(aa27::A, aa27::C));
    EXPECT_EQ(-3,    scheme.score(aa27::C, aa27::A));
    EXPECT_EQ( 9,    scheme.score(aa27::D, aa27::D));
    EXPECT_EQ( 0,    scheme.score(aa27::N, aa27::A));

    // Test set function
    scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM45);

    EXPECT_EQ( 5,    scheme.score(aa27::A, aa27::A));
    EXPECT_EQ(-1,    scheme.score(aa27::A, aa27::C));
    EXPECT_EQ(-1,    scheme.score(aa27::C, aa27::A));
    EXPECT_EQ( 7,    scheme.score(aa27::D, aa27::D));
    EXPECT_EQ(-1,    scheme.score(aa27::N, aa27::A));

    scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM62);

    EXPECT_EQ( 4,    scheme.score(aa27::A, aa27::A));
    EXPECT_EQ( 0,    scheme.score(aa27::A, aa27::C));
    EXPECT_EQ( 0,    scheme.score(aa27::C, aa27::A));
    EXPECT_EQ( 6,    scheme.score(aa27::D, aa27::D));
    EXPECT_EQ(-2,    scheme.score(aa27::N, aa27::A));

    scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM80);

    EXPECT_EQ( 7,    scheme.score(aa27::A, aa27::A));
    EXPECT_EQ(-1,    scheme.score(aa27::A, aa27::C));
    EXPECT_EQ(-1,    scheme.score(aa27::C, aa27::A));
    EXPECT_EQ(10,    scheme.score(aa27::D, aa27::D));
    EXPECT_EQ(-3,    scheme.score(aa27::N, aa27::A));
}
