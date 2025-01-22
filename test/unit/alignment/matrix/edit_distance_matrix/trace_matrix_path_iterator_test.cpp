// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/edit_distance_trace_matrix_full.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_concept.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "edit_distance_trace_matrix.hpp"

struct trace_iterator_fixture : public ::testing::Test
{
    using sequence_t = std::string;
    using trace_path_vector = std::vector<seqan3::detail::trace_directions>;

    void SetUp() override
    {
        matrix.reserve(10u);

        matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                          {0b0000'0000u, 0b0000'0000u, 0b0u},
                          {0b1111'1111u, 0b1111'1111u, 0b1u});

        matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                          {0b1000'0011u, 0b0000'0011u, 0b1u},
                          {0b1111'1110u, 0b1111'1111u, 0b1u});

        matrix.add_column({0b0000'0001u, 0b0000'0000u, 0b0u},
                          {0b0000'1110u, 0b0000'1100u, 0b0u},
                          {0b1111'1000u, 0b1111'1111u, 0b1u});

        matrix.add_column({0b0000'0111u, 0b0000'0000u, 0b0u},
                          {0b0011'1110u, 0b0011'0000u, 0b0u},
                          {0b1110'0000u, 0b1111'1111u, 0b1u});

        matrix.add_column({0b0001'1111u, 0b0000'0000u, 0b0u},
                          {0b0111'1110u, 0b1100'0000u, 0b0u},
                          {0b1000'0000u, 0b1111'1111u, 0b1u});

        matrix.add_column({0b0111'1101u, 0b0000'0000u, 0b0u},
                          {0b1111'1111u, 0b0000'0011u, 0b1u},
                          {0b0000'0100u, 0b1111'1111u, 0b1u});

        matrix.add_column({0b1111'0011u, 0b0000'0000u, 0b0u},
                          {0b0111'1100u, 0b0000'1111u, 0b0u},
                          {0b0001'1000u, 0b1111'1010u, 0b1u});

        matrix.add_column({0b1100'0111u, 0b0000'0101u, 0b0u},
                          {0b0111'1000u, 0b0011'1111u, 0b0u},
                          {0b0110'0000u, 0b1110'0100u, 0b1u});

        matrix.add_column({0b1001'1111u, 0b0001'1011u, 0b0u},
                          {0b0111'1000u, 0b1111'1111u, 0b0u},
                          {0b1000'0000u, 0b1000'1000u, 0b1u});

        matrix.add_column({0b0111'1111u, 0b0111'0100u, 0b0u},
                          {0b1111'1011u, 0b1111'1111u, 0b1u},
                          {0b0000'0000u, 0b0001'0101u, 0b0u});
    }

    using expect_matrix_type = seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions>;
    // This is nearly seqan3::test::alignment::fixture::global::edit_distance::unbanded::dna4_02T, but the second
    // sequence has an additional A.
    expect_matrix_type expect_matrix{seqan3::detail::number_rows{18u},
                                     seqan3::detail::number_cols{10u},
                                     {//     e,  A,  C,  G,  T,  A,  C,  G,  T,  A
                                      /*e*/ N, l,  l,  l,  l,  l,   l,   l,   l,   l,
                                      /*A*/ u, D,  l,  l,  l,  Dl,  l,   l,   l,   Dl,
                                      /*A*/ u, Du, D,  Dl, Dl, D,   l,   l,   l,   Dl,
                                      /*C*/ u, u,  D,  Dl, Dl, Dul, D,   l,   l,   l,
                                      /*C*/ u, u,  Du, D,  Dl, Dl,  Du,  D,   Dl,  Dl,
                                      /*G*/ u, u,  u,  D,  Dl, Dl,  Dul, D,   Dl,  Dl,
                                      /*G*/ u, u,  u,  Du, D,  Dl,  Dl,  Du,  D,   Dl,
                                      /*T*/ u, u,  u,  u,  D,  Dl,  Dl,  Dul, D,   Dl,
                                      /*T*/ u, Du, u,  u,  u,  D,   l,   l,   ul,  D,
                                      /*A*/ u, Du, u,  u,  u,  Du,  D,   Dl,  Dl,  Du,
                                      /*A*/ u, Du, u,  u,  u,  Du,  Du,  D,   Dl,  D,
                                      /*C*/ u, u,  Du, u,  u,  u,   D,   Dul, D,   Dul,
                                      /*C*/ u, u,  Du, u,  u,  u,   Du,  D,   Dul, D,
                                      /*G*/ u, u,  u,  Du, u,  u,   u,   D,   Dl,  Dul,
                                      /*G*/ u, u,  u,  Du, u,  u,   u,   Du,  D,   Dl,
                                      /*T*/ u, u,  u,  u,  Du, u,   u,   u,   D,   Dl,
                                      /*T*/ u, u,  u,  u,  Du, u,   u,   u,   Du,  D,
                                      /*A*/ u, Du, u,  u,  u,  Du,  u,   u,   u,   D}};

    auto path(size_t row, size_t column)
    {
        return matrix.trace_path({seqan3::detail::row_index_type{row}, seqan3::detail::column_index_type{column}});
    };

    sequence_t sequence1 = "ACGTACGTA";
    sequence_t sequence2 = "AACCGGTAAACCGGTTA";
    seqan3::detail::aligned_sequence_builder<sequence_t &, sequence_t &> builder{sequence1, sequence2};

    matrix_type<false, false> matrix{18u};
};

TEST_F(trace_iterator_fixture, trace_matrix)
{
    EXPECT_EQ(matrix, expect_matrix);
}

TEST_F(trace_iterator_fixture, trace_paths)
{
    EXPECT_RANGE_EQ(path(0u, 0u), (trace_path_vector{}));
    EXPECT_RANGE_EQ(path(1u, 1u), (trace_path_vector{D}));
    EXPECT_RANGE_EQ(path(0u, 9u), (trace_path_vector{l, l, l, l, l, l, l, l, l}));
    EXPECT_RANGE_EQ(path(17u, 0u), (trace_path_vector{u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u}));
    EXPECT_RANGE_EQ(path(1u, 2u), (trace_path_vector{l, D}));
    EXPECT_RANGE_EQ(path(2u, 1u), (trace_path_vector{/*D*/ u, D}));
    EXPECT_RANGE_EQ(path(7u, 9u), (trace_path_vector{/*D*/ l, D, /*D*/ u, D, /*D*/ u, D, D, l, l, l, D}));
    EXPECT_RANGE_EQ(path(10u, 9u),
                    (trace_path_vector{D, /*D*/ l, /*D*/ l, D, D, D, /*D*/ u, D, /*D*/ u, D, /*D*/ u, D}));
    EXPECT_RANGE_EQ(path(11u, 9u), (trace_path_vector{/*Du*/ l, D, D, D, D, D, /*D*/ u, D, /*D*/ u, D, /*D*/ u, D}));
    EXPECT_RANGE_EQ(path(17u, 9u),
                    (trace_path_vector{D,
                                       /*D*/ u,
                                       D,
                                       /*D*/ u,
                                       D,
                                       /*D*/ u,
                                       D,
                                       /*D*/ u,
                                       /*D*/ u,
                                       D,
                                       D,
                                       /*D*/ u,
                                       D,
                                       /*D*/ u,
                                       D,
                                       /*D*/ u,
                                       D}));
}

TEST_F(trace_iterator_fixture, invalid_trace_path)
{
    EXPECT_EQ(matrix.cols(), 10u);
    EXPECT_EQ(matrix.rows(), 18u);

    EXPECT_THROW(path(0u, 10u), std::invalid_argument);
    EXPECT_THROW(path(18u, 0u), std::invalid_argument);
    EXPECT_THROW(path(18u, 9u), std::invalid_argument);
    EXPECT_THROW(path(17u, 10u), std::invalid_argument);
}

TEST_F(trace_iterator_fixture, alignment_0_0)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(0u, 0u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, ""s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, ""s);
}

TEST_F(trace_iterator_fixture, alignment_1_1)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(1u, 1u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "A"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "A"s);
}

TEST_F(trace_iterator_fixture, alignment_0_9)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(0u, 9u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "ACGTACGTA"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "---------"s);
}

TEST_F(trace_iterator_fixture, alignment_17_0)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(17u, 0u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "-----------------"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "AACCGGTAAACCGGTTA"s);
}

TEST_F(trace_iterator_fixture, alignment_1_2)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(1u, 2u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "AC"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "A-"s);
}

TEST_F(trace_iterator_fixture, alignment_2_1)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(2u, 1u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "A-"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "AA"s);
}

TEST_F(trace_iterator_fixture, alignment_7_9)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(7u, 9u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "ACGTAC-G-TA"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "A---ACCGGT-"s);
}

TEST_F(trace_iterator_fixture, alignment_10_9)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(10u, 9u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "A-C-G-TACGTA"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "AACCGGTAA--A"s);
}

TEST_F(trace_iterator_fixture, alignment_11_9)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(11u, 9u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "A-C-G-TACGTA"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "AACCGGTAAAC-"s);
}

TEST_F(trace_iterator_fixture, alignment_17_9)
{
    using namespace std::string_literals;

    auto alignment_result = builder(path(17u, 9u));
    EXPECT_RANGE_EQ(std::get<0>(alignment_result.alignment) | seqan3::views::to_char, "A-C-G-TA--C-G-T-A"s);
    EXPECT_RANGE_EQ(std::get<1>(alignment_result.alignment) | seqan3::views::to_char, "AACCGGTAAACCGGTTA"s);
}
