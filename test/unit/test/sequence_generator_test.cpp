// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

// generate (single) fixed sized sequence
TEST(random_sequence_generator, fixed_length)
{
    std::mt19937_64 random_engine{0};

    seqan3::test::random_sequence_generator<seqan3::dna4_vector> random_sequence_generator{3};

#ifdef _LIBCPP_VERSION
    EXPECT_EQ(random_sequence_generator(random_engine), "GTC"_dna4);
    EXPECT_EQ(random_sequence_generator(random_engine), "GAG"_dna4);
#else
    EXPECT_EQ(random_sequence_generator(random_engine), "TAG"_dna4);
    EXPECT_EQ(random_sequence_generator(random_engine), "AGC"_dna4);
#endif
}

// generate (single) variable sized sequence, generates sequence with size ± size_variance.
TEST(random_sequence_generator, variable_length_and_different_random_engines)
{
    // Test different uniform random bit generators
    std::ranlux24 random_engine1{1};
    std::ranlux48 random_engine2{1};
    std::mt19937_64 random_engine3{1};

    seqan3::test::random_sequence_generator<seqan3::dna5_vector> random_sequence_generator{4, 2};

#ifdef _LIBCPP_VERSION
    EXPECT_EQ(random_sequence_generator(random_engine1), "TCCGTT"_dna5);
    EXPECT_EQ(random_sequence_generator(random_engine2), "CGTAAG"_dna5);
    EXPECT_EQ(random_sequence_generator(random_engine3), "GA"_dna5);
#else
    EXPECT_EQ(random_sequence_generator(random_engine1), "CCAN"_dna5);
    EXPECT_EQ(random_sequence_generator(random_engine2), "NN"_dna5);
    EXPECT_EQ(random_sequence_generator(random_engine3), "AG"_dna5);
#endif
}

// generate a collection of variable/fixed sized sequences
TEST(random_sequence_generator, sequence_collection)
{
    seqan3::test::random_sequence_generator<seqan3::dna5_vector> random_sequence_generator{4, 2};
    std::mt19937_64 random_engine{0};

    std::vector<seqan3::dna5_vector> sequences(4);

    std::ranges::generate(sequences,
                          [&]()
                          {
                              return random_sequence_generator(random_engine);
                          });

#ifdef _LIBCPP_VERSION
    std::vector<seqan3::dna5_vector> resulting_sequences = {"CTACN"_dna5, "AC"_dna5, "CTCGG"_dna5, "GGCANN"_dna5};
#else
    std::vector<seqan3::dna5_vector> resulting_sequences = {"TA"_dna5, "GANG"_dna5, "TGNTCC"_dna5, "NNGC"_dna5};
#endif
    EXPECT_EQ(sequences, resulting_sequences);
}

// generate sequence pairs
TEST(random_sequence_generator, sequence_pairs)
{
    seqan3::test::random_sequence_generator<seqan3::dna5_vector> random_sequence_generator{4, 2};
    std::mt19937_64 random_engine{0};

    std::vector<std::tuple<seqan3::dna5_vector, seqan3::dna5_vector>> sequence_pairs(3);

    std::ranges::generate(sequence_pairs,
                          [&]() -> std::tuple<seqan3::dna5_vector, seqan3::dna5_vector>
                          {
                              return {random_sequence_generator(random_engine),
                                      random_sequence_generator(random_engine)};
                          });

#ifdef _LIBCPP_VERSION
    std::vector<std::tuple<seqan3::dna5_vector, seqan3::dna5_vector>> resulting_pairs{{"CTACN"_dna5, "AC"_dna5},
                                                                                      {"CTCGG"_dna5, "GGCANN"_dna5},
                                                                                      {"CACA"_dna5, "CNC"_dna5}};
#else
    std::vector<std::tuple<seqan3::dna5_vector, seqan3::dna5_vector>> resulting_pairs{{"TA"_dna5, "GANG"_dna5},
                                                                                      {"TGNTCC"_dna5, "NNGC"_dna5},
                                                                                      {"GG"_dna5, "AA"_dna5}};
#endif

    EXPECT_EQ(sequence_pairs, resulting_pairs);
}
