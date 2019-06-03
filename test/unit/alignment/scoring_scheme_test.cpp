// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
#include <seqan3/test/cereal.hpp>

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
    EXPECT_TRUE((ScoringScheme<TypeParam, alph_t>));
    EXPECT_TRUE((ScoringScheme<TypeParam const, alph_t>));
    EXPECT_TRUE((ScoringScheme<TypeParam const &, alph_t>));
    EXPECT_FALSE((ScoringScheme<TypeParam const &, char>));
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

    EXPECT_EQ(matrix_size, alphabet_size<alph_t>);
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

    for (uint8_t i = 0; i < alphabet_size<alph_t>; ++i)
    {
        for (uint8_t j = 0; j < alphabet_size<alph_t>; ++j)
        {
            int8_t expected = i == j ? 5 : -3;
            EXPECT_EQ(expected, scheme.score(assign_rank_to(i, alph_t{}), assign_rank_to(j, alph_t{})));
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

    for (uint8_t i = 0; i < alphabet_size<alph_t>; ++i)
    {
        for (uint8_t j = 0; j < alphabet_size<alph_t>; ++j)
        {
            int8_t expected = i == j ? 0 : -1;
            EXPECT_EQ(expected, scheme.score(assign_rank_to(i, alph_t{}), assign_rank_to(j, alph_t{})));
        }
    }
}

TYPED_TEST(generic, custom)
{
    using alph_t = typename TestFixture::alph_t;

    typename TypeParam::matrix_type matrix;

    for (uint8_t i = 0; i < alphabet_size<alph_t>; ++i)
        for (uint8_t j = 0; j < alphabet_size<alph_t>; ++j)
            matrix[i][j] = i * i + j;

    // Test constructor
    TypeParam scheme{matrix};
    // Test set function
    scheme.set_custom_matrix(matrix);

    if constexpr (detail::is_type_specialisation_of_v<TypeParam, aminoacid_scoring_scheme>)
    {
        EXPECT_EQ(0*0+0,    scheme.score('A'_aa27, 'A'_aa27));
        EXPECT_EQ(0*0+2,    scheme.score('A'_aa27, 'C'_aa27));
        EXPECT_EQ(2*2+0,    scheme.score('C'_aa27, 'A'_aa27));
        EXPECT_EQ(8*8+8,    scheme.score('I'_aa27, 'I'_aa27));
        EXPECT_EQ(0*0+13,   scheme.score('A'_aa27, 'N'_aa27));
        EXPECT_EQ(2*2+1,    scheme.score('C'_aa27, 'B'_aa27));
    } else
    {
        EXPECT_EQ(0*0+0,    scheme.score('A'_dna15, 'A'_dna15)); // A is 0th
        EXPECT_EQ(0*0+2,    scheme.score('A'_dna15, 'C'_dna15)); // C is 2nd
        EXPECT_EQ(2*2+0,    scheme.score('C'_dna15, 'A'_dna15));
        EXPECT_EQ(3*3+3,    scheme.score('D'_dna15, 'D'_dna15)); // D is 3rd
        EXPECT_EQ(8*8+0,    scheme.score('N'_dna15, 'A'_dna15)); // N is 8th
        EXPECT_EQ(0*0+8,    scheme.score('A'_dna15, 'N'_dna15));
    }
}

TYPED_TEST(generic, convertability)
{
    using alph_t = typename TestFixture::alph_t;

    typename TypeParam::matrix_type matrix;

    for (uint8_t i = 0; i < alphabet_size<alph_t>; ++i)
        for (uint8_t j = 0; j < alphabet_size<alph_t>; ++j)
            matrix[i][j] = i * i + j;

    TypeParam scheme{};
    scheme.set_custom_matrix(matrix);

    if constexpr (detail::is_type_specialisation_of_v<TypeParam, aminoacid_scoring_scheme>)
    {
        using aa_types = type_list<aa27, aa20>;
        meta::for_each(aa_types{}, [&] (auto && aa) constexpr
        {
            using nucl_t = std::decay_t<decltype(aa)>;

            EXPECT_EQ(scheme.score('C'_aa27,                  'G'_aa27),
                      scheme.score(nucl_t{}.assign_char('C'), nucl_t{}.assign_char('G')));
            EXPECT_EQ(scheme.score('T'_aa27,                  nucl_t{}.assign_char('T')),
                      scheme.score(nucl_t{}.assign_char('T'), 'T'_aa27));
        });
    } else
    {

        using nucl_types = type_list<dna4, dna5, dna15, rna4, rna5, rna15>;
        meta::for_each(nucl_types{}, [&] (auto && nucl) constexpr
        {
            using nucl_t = std::decay_t<decltype(nucl)>;

            EXPECT_EQ(scheme.score('C'_dna15,                 'G'_dna15),
                      scheme.score(nucl_t{}.assign_char('C'), nucl_t{}.assign_char('G')));
            EXPECT_EQ(scheme.score('T'_dna15,                 'T'_dna15),
                      scheme.score(nucl_t{}.assign_char('T'), nucl_t{}.assign_char('T')));
            EXPECT_EQ(scheme.score('A'_dna15,                 'C'_dna15),
                      scheme.score(nucl_t{}.assign_char('A'), nucl_t{}.assign_char('C')));

            EXPECT_EQ(scheme.score('C'_dna15,                 nucl_t{}.assign_char('G')),
                      scheme.score(nucl_t{}.assign_char('C'), 'G'_dna15));
            EXPECT_EQ(scheme.score('C'_dna15,                 nucl_t{}.assign_char('A')),
                      scheme.score(nucl_t{}.assign_char('C'), 'A'_dna15));
        });

    }
}

TYPED_TEST(generic, serialisation)
{
    TypeParam scheme1;

    scheme1.set_hamming_distance();
    test::do_serialisation(scheme1);

    scheme1.set_simple_scheme(match_score{11}, mismatch_score{-7});
    test::do_serialisation(scheme1);
}

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
    EXPECT_EQ( 4,    scheme.score('A'_aa27, 'A'_aa27));
    EXPECT_EQ(-3,    scheme.score('A'_aa27, 'C'_aa27));
    EXPECT_EQ(-3,    scheme.score('C'_aa27, 'A'_aa27));
    EXPECT_EQ( 9,    scheme.score('D'_aa27, 'D'_aa27));
    EXPECT_EQ( 0,    scheme.score('N'_aa27, 'A'_aa27));

    // Test set function
    scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM45);

    EXPECT_EQ( 5,    scheme.score('A'_aa27, 'A'_aa27));
    EXPECT_EQ(-1,    scheme.score('A'_aa27, 'C'_aa27));
    EXPECT_EQ(-1,    scheme.score('C'_aa27, 'A'_aa27));
    EXPECT_EQ( 7,    scheme.score('D'_aa27, 'D'_aa27));
    EXPECT_EQ(-1,    scheme.score('N'_aa27, 'A'_aa27));

    scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM62);

    EXPECT_EQ( 4,    scheme.score('A'_aa27, 'A'_aa27));
    EXPECT_EQ( 0,    scheme.score('A'_aa27, 'C'_aa27));
    EXPECT_EQ( 0,    scheme.score('C'_aa27, 'A'_aa27));
    EXPECT_EQ( 6,    scheme.score('D'_aa27, 'D'_aa27));
    EXPECT_EQ(-2,    scheme.score('N'_aa27, 'A'_aa27));

    scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM80);

    EXPECT_EQ( 7,    scheme.score('A'_aa27, 'A'_aa27));
    EXPECT_EQ(-1,    scheme.score('A'_aa27, 'C'_aa27));
    EXPECT_EQ(-1,    scheme.score('C'_aa27, 'A'_aa27));
    EXPECT_EQ(10,    scheme.score('D'_aa27, 'D'_aa27));
    EXPECT_EQ(-3,    scheme.score('N'_aa27, 'A'_aa27));
}
