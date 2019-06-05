// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/edit_distance_trace_matrix_full.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "../pairwise/fixture/alignment_fixture.hpp"

using namespace seqan3;
using namespace seqan3::test::alignment::fixture;

using word_type = uint8_t;

template <bool is_semi_global, bool use_max_errors>
class matrix_type
    : public seqan3::detail::edit_distance_trace_matrix_full<word_type, is_semi_global, use_max_errors>
{
public:
    using base_t = seqan3::detail::edit_distance_trace_matrix_full<word_type, is_semi_global, use_max_errors>;

    matrix_type(size_t const rows_size) : base_t{rows_size}
    {}

    using base_t::add_column;
    using base_t::reserve;
};

std::vector<std::vector<detail::trace_directions>> as_row_wise_vector(auto matrix)
{
    std::vector<std::vector<detail::trace_directions>> result{};
    for (unsigned row = 0; row < matrix.rows(); ++row)
    {
        result.push_back({});
        for (unsigned col = 0; col < matrix.cols(); ++col)
            result.back().push_back(matrix.at(row, col));
    }
    return result;
}

TEST(global, empty)
{
    matrix_type<false, false> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(global, epsilon)
{
    matrix_type<false, false> matrix{1u};

    matrix.add_column({}, {}, {});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON}};

    EXPECT_EQ(result, expect);
}

TEST(global, epsilon_row)
{
    matrix_type<false, false> matrix{1u};

    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON, L, L, L, L}};

    EXPECT_EQ(result, expect);
}

TEST(global, single_word)
{
    matrix_type<false, false> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u});
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u});
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u});
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u});
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u});
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u});
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u});
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u});
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u});
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  },
        {U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL },
        {U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  },
        {U  ,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  },
        {U  ,U  ,U  ,D  ,DL ,DUL,DU ,DU ,DU ,DU },
        {U  ,U  ,U  ,U  ,D  ,D  ,DL ,DUL,DU ,DU },
        {U  ,U  ,U  ,U  ,DU ,DU ,D  ,D  ,DL ,DUL}
    };

    EXPECT_EQ(result, expect);
}

TEST(global, multiple_words)
{
    matrix_type<false, false> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0011u, 0b0000'0011u, 0b0u},
                      {0b1111'1110u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0001u, 0b0000'0000u, 0b0u},
                      {0b0000'1110u, 0b0000'1100u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0111u, 0b0000'0000u, 0b0u},
                      {0b0011'1110u, 0b0011'0000u, 0b0u},
                      {0b1110'0000u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0001'1111u, 0b0000'0000u, 0b0u},
                      {0b1111'1110u, 0b1100'0000u, 0b1u},
                      {0b1000'0000u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0111'1101u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b0000'0011u, 0b0u},
                      {0b0000'0100u, 0b1111'1110u, 0b1u});
    matrix.add_column({0b1111'0011u, 0b0000'0001u, 0b0u},
                      {0b1111'1100u, 0b0000'1111u, 0b0u},
                      {0b0001'1000u, 0b1111'1000u, 0b1u});
    matrix.add_column({0b1100'0111u, 0b0000'0111u, 0b0u},
                      {0b1111'1000u, 0b0011'1111u, 0b0u},
                      {0b0110'0000u, 0b1110'0000u, 0b1u});
    matrix.add_column({0b0001'1111u, 0b0001'1111u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u},
                      {0b1000'0000u, 0b1000'0001u, 0b1u});
    matrix.add_column({0b0111'1111u, 0b0111'1100u, 0b0u},
                      {0b1111'1011u, 0b1111'1111u, 0b1u},
                      {0b0000'0000u, 0b0000'0110u, 0b0u});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  },
        {U  ,D  ,L  ,L  ,L  ,DL ,L  ,L  ,L  ,DL },
        {U  ,DU ,D  ,DL ,DL ,D  ,L  ,L  ,L  ,DL },
        {U  ,U  ,D  ,DL ,DL ,DUL,D  ,L  ,L  ,L  },
        {U  ,U  ,DU ,D  ,DL ,DL ,DU ,D  ,DL ,DL },
        {U  ,U  ,U  ,D  ,DL ,DL ,DUL,D  ,DL ,DL },
        {U  ,U  ,U  ,DU ,D  ,DL ,DL ,DU ,D  ,DL },
        {U  ,U  ,U  ,U  ,D  ,DL ,DL ,DUL,D  ,DL },
        {U  ,U  ,U  ,U  ,DU ,D  ,DL ,DL ,DU ,D  },
        {U  ,DU ,U  ,U  ,U  ,D  ,DL ,DL ,DUL,D  },
        {U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,DL ,DU },
        {U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DL ,DUL},
        {U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,DL },
        {U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DL },
        {U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL },
        {U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL },
        {U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  },
        {U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  }
    };

    EXPECT_EQ(result, expect);
}

TEST(semi_global, empty)
{
    matrix_type<true, false> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global, epsilon)
{
    matrix_type<true, false> matrix{1u};

    matrix.add_column({}, {}, {});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON}};

    EXPECT_EQ(result, expect);
}
//
TEST(semi_global, epsilon_row)
{
    matrix_type<true, false> matrix{1u};

    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON, NON, NON, NON, NON}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global, single_word)
{
    matrix_type<true, false> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u});
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u});
    matrix.add_column({0b0000'0000u}, {0b0001'1111u}, {0b1110'1110u});
    matrix.add_column({0b0000'0001u}, {0b0010'0011u}, {0b1101'1101u});
    matrix.add_column({0b0000'0010u}, {0b1111'1111u}, {0b1101'1001u});
    matrix.add_column({0b0010'0010u}, {0b0111'1111u}, {0b1011'1011u});
    matrix.add_column({0b0100'0100u}, {0b1111'1111u}, {0b0011'0011u});
    matrix.add_column({0b0100'0100u}, {0b1111'1111u}, {0b0111'0111u});
    matrix.add_column({0b1000'1000u}, {0b1111'1111u}, {0b0110'0111u});
    matrix.add_column({0b1000'0000u}, {0b1111'0001u}, {0b1110'1110u});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON},
        {U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,D  },
        {U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,U  },
        {U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,U  },
        {U  ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,U  },
        {U  ,DU ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,D  },
        {U  ,U  ,U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU },
        {U  ,U  ,U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU },
        {U  ,U  ,U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL}
    };

    EXPECT_EQ(result, expect);
}

TEST(semi_global, multiple_words)
{
    matrix_type<true, false> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0011u, 0b0000'0011u, 0b0u},
                      {0b1111'1110u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0001u, 0b0000'0000u, 0b0u},
                      {0b0000'1111u, 0b0000'1100u, 0b0u},
                      {0b1111'1001u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0110u, 0b0000'0000u, 0b0u},
                      {0b0011'1111u, 0b0011'0000u, 0b0u},
                      {0b1110'0011u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0001'1100u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1100'0000u, 0b1u},
                      {0b1000'0111u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0111'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'0011u, 0b0000'0011u, 0b0u},
                      {0b0001'1110u, 0b1111'1110u, 0b1u});
    matrix.add_column({0b1100'0001u, 0b0000'0001u, 0b0u},
                      {0b1100'1111u, 0b0000'1111u, 0b0u},
                      {0b0111'1101u, 0b1111'1000u, 0b1u});
    matrix.add_column({0b0000'0010u, 0b0000'0111u, 0b0u},
                      {0b0011'1111u, 0b0011'1111u, 0b0u},
                      {0b1111'0001u, 0b1110'0001u, 0b1u});
    matrix.add_column({0b0000'1110u, 0b0001'1100u, 0b0u},
                      {0b1111'1111u, 0b1111'1100u, 0b1u},
                      {0b1100'0011u, 0b1000'0111u, 0b1u});
    matrix.add_column({0b0000'1000u, 0b0111'0000u, 0b0u},
                      {0b1111'0011u, 0b1111'0011u, 0b1u},
                      {0b0100'1110u, 0b0001'1111u, 0b0u});

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON},
        {U  ,D  ,DUL,DU ,DU ,D  ,DUL,DU ,DU ,D  },
        {U  ,DU ,D  ,DUL,DU ,DU ,D  ,DL ,DUL,DU },
        {U  ,U  ,D  ,DL ,DUL,U  ,DU ,D  ,DL ,U  },
        {U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,DL ,UL },
        {U  ,U  ,U  ,D  ,DL ,DUL,U  ,DU ,D  ,D  },
        {U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,D  },
        {U  ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,DU ,DU },
        {U  ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  },
        {U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,DU },
        {U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU },
        {U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,U  },
        {U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,U  },
        {U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL},
        {U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL },
        {U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL },
        {U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  },
        {U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  }
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, empty)
{
    matrix_type<false, true> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, epsilon)
{
    matrix_type<false, true> matrix{1u};

    matrix.add_column({}, {}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, epsilon_row)
{
    matrix_type<false, true> matrix{1u};

    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 0u);
    matrix.add_column({}, {}, {}, 0u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON, L, L,NON,NON}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_1)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u}, 6u);
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u}, 7u);
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u}, 8u);
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u}, 9u);
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u}, 9u);
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u}, 9u);
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u}, 9u);
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u}, 9u);
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u}, 7u);
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u}, 7u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  },
        {U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL },
        {U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  },
        {U  ,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  },
        {NON,U  ,U  ,D  ,DL ,DUL,DU ,DU ,DU ,DU },
        {NON,NON,U  ,U  ,D  ,D  ,DL ,DUL,NON,NON},
        {NON,NON,NON,U  ,DU ,DU ,D  ,D  ,NON,NON}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_2)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u}, 5u);
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u}, 6u);
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u}, 7u);
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u}, 8u);
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u}, 8u);
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u}, 8u);
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u}, 6u);
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u}, 6u);
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u}, 6u);
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u}, 6u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  },
        {U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL },
        {U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  },
        {NON,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  },
        {NON,NON,U  ,D  ,DL ,DUL,NON,NON,NON,NON},
        {NON,NON,NON,U  ,D  ,D  ,NON,NON,NON,NON},
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_3)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u}, 4u);
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u}, 5u);
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u}, 6u);
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u}, 7u);
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u}, 5u);
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u}, 5u);
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u}, 5u);
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u}, 5u);
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u}, 0u);
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u}, 0u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,NON,NON},
        {U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,NON,NON},
        {U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,NON,NON},
        {U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,NON,NON},
        {NON,U  ,DU ,DU ,DU ,DU ,D  ,D  ,NON,NON},
        {NON,NON,D  ,DUL,NON,NON,NON,NON,NON,NON},
        {NON,NON,NON,D  ,NON,NON,NON,NON,NON,NON},
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON},
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, multiple_words_1)
{
    matrix_type<false, true> matrix{10u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0u}, {0b0000'0000u, 0b0u}, {0b1111'1111u, 0b1u}, 6u);
    matrix.add_column({0b0000'0000u, 0b0u}, {0b0001'0001u, 0b1u}, {0b1111'1110u, 0b1u}, 7u);
    matrix.add_column({0b0000'0001u, 0b0u}, {0b0001'1111u, 0b1u}, {0b1110'1100u, 0b1u}, 8u);
    matrix.add_column({0b0001'0001u, 0b0u}, {0b0011'1110u, 0b0u}, {0b1101'1100u, 0b1u}, 9u);
    matrix.add_column({0b0010'0011u, 0b0u}, {0b1111'1110u, 0b1u}, {0b1001'1000u, 0b1u}, 9u);
    matrix.add_column({0b0010'0011u, 0b0u}, {0b1111'1100u, 0b1u}, {0b1011'1000u, 0b1u}, 9u);
    matrix.add_column({0b0100'0111u, 0b0u}, {0b1111'1100u, 0b1u}, {0b0011'0000u, 0b1u}, 9u);
    matrix.add_column({0b0100'0111u, 0b0u}, {0b1111'1000u, 0b1u}, {0b0111'0000u, 0b1u}, 9u);
    matrix.add_column({0b1000'1111u, 0b0u}, {0b1111'1000u, 0b1u}, {0b0110'0000u, 0b0u}, 7u);
    matrix.add_column({0b1000'1111u, 0b0u}, {0b1111'0001u, 0b1u}, {0b1110'0000u, 0b0u}, 7u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  },
        {U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL },
        {U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  },
        {U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  },
        {U  ,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  },
        {NON,U  ,U  ,D  ,DL ,DUL,DU ,DU ,DU ,DU },
        {NON,NON,U  ,U  ,D  ,D  ,DL ,DUL,NON,NON},
        {NON,NON,NON,U  ,DU ,DU ,D  ,D  ,NON,NON},
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, multiple_words_2)
{
    matrix_type<false, true> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1111'1111u, 0b1u}, 9u);
    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0011u, 0b0000'0011u, 0b0u},
                      {0b1111'1110u, 0b1111'1111u, 0b1u}, 10u);
    matrix.add_column({0b0000'0001u, 0b0000'0000u, 0b0u},
                      {0b0000'1110u, 0b0000'1100u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u}, 11u);
    matrix.add_column({0b0000'0111u, 0b0000'0000u, 0b0u},
                      {0b0011'1110u, 0b0011'0000u, 0b0u},
                      {0b1110'0000u, 0b1111'1111u, 0b1u}, 12u);
    matrix.add_column({0b0001'1111u, 0b0000'0000u, 0b0u},
                      {0b1111'1110u, 0b1100'0000u, 0b1u},
                      {0b1000'0000u, 0b1111'1111u, 0b1u}, 13u);
    matrix.add_column({0b0111'1101u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b0000'0011u, 0b0u},
                      {0b0000'0100u, 0b1111'1110u, 0b1u}, 14u);
    matrix.add_column({0b1111'0011u, 0b0000'0001u, 0b0u},
                      {0b1111'1100u, 0b0000'1111u, 0b0u},
                      {0b0001'1000u, 0b1111'1000u, 0b1u}, 15u);
    matrix.add_column({0b1100'0111u, 0b0000'0111u, 0b0u},
                      {0b1111'1000u, 0b0011'1111u, 0b0u},
                      {0b0110'0000u, 0b1110'0000u, 0b1u}, 16u);
    matrix.add_column({0b0001'1111u, 0b0001'1111u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u},
                      {0b1000'0000u, 0b1000'0001u, 0b1u}, 17u);
    matrix.add_column({0b0111'1111u, 0b0111'1100u, 0b0u},
                      {0b1111'1011u, 0b1111'1111u, 0b1u},
                      {0b0000'0000u, 0b0000'0110u, 0b0u}, 18u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  },
        {U  ,D  ,L  ,L  ,L  ,DL ,L  ,L  ,L  ,DL },
        {U  ,DU ,D  ,DL ,DL ,D  ,L  ,L  ,L  ,DL },
        {U  ,U  ,D  ,DL ,DL ,DUL,D  ,L  ,L  ,L  },
        {U  ,U  ,DU ,D  ,DL ,DL ,DU ,D  ,DL ,DL },
        {U  ,U  ,U  ,D  ,DL ,DL ,DUL,D  ,DL ,DL },
        {U  ,U  ,U  ,DU ,D  ,DL ,DL ,DU ,D  ,DL },
        {U  ,U  ,U  ,U  ,D  ,DL ,DL ,DUL,D  ,DL },
        {U  ,U  ,U  ,U  ,DU ,D  ,DL ,DL ,DU ,D  },
        {NON,DU ,U  ,U  ,U  ,D  ,DL ,DL ,DUL,D  },
        {NON,NON,U  ,U  ,U  ,DU ,D  ,DL ,DL ,DU },
        {NON,NON,NON,U  ,U  ,U  ,D  ,DL ,DL ,DUL},
        {NON,NON,NON,NON,U  ,U  ,DU ,D  ,DL ,DL },
        {NON,NON,NON,NON,NON,U  ,U  ,D  ,DL ,DL },
        {NON,NON,NON,NON,NON,NON,U  ,DU ,D  ,DL },
        {NON,NON,NON,NON,NON,NON,NON,U  ,D  ,DL },
        {NON,NON,NON,NON,NON,NON,NON,NON,DU ,D  },
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,D  }
    };

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, empty)
{
    matrix_type<true, true> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, epsilon)
{
    matrix_type<true, true> matrix{1u};

    matrix.add_column({}, {}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, epsilon_row)
{
    matrix_type<true, true> matrix{1u};

    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect{{NON,NON,NON,NON,NON}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, single_word)
{
    matrix_type<true, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u}, 6u);
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u}, 7u);
    matrix.add_column({0b0000'0000u}, {0b0001'1111u}, {0b1110'1110u}, 8u);
    matrix.add_column({0b0000'0001u}, {0b0010'0011u}, {0b1101'1101u}, 9u);
    matrix.add_column({0b0000'0010u}, {0b1111'1111u}, {0b1101'1001u}, 9u);
    matrix.add_column({0b0010'0010u}, {0b0111'1111u}, {0b1011'1011u}, 9u);
    matrix.add_column({0b0100'0100u}, {0b1111'1111u}, {0b0011'0011u}, 9u);
    matrix.add_column({0b0100'0100u}, {0b1111'1111u}, {0b0111'0111u}, 9u);
    matrix.add_column({0b1000'1000u}, {0b1111'1111u}, {0b0110'0111u}, 9u);
    matrix.add_column({0b1000'0000u}, {0b1111'0001u}, {0b1110'1110u}, 8u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON},
        {U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,D  },
        {U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,U  },
        {U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,U  },
        {U  ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,U  },
        {U  ,DU ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,D  },
        {NON,U  ,U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU },
        {NON,NON,U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU },
        {NON,NON,NON,U  ,DU ,U  ,D  ,D  ,DL ,NON}
    };

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, multiple_words)
{
    matrix_type<true, true> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1111'1111u, 0b1u}, 9u);
    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0011u, 0b0000'0011u, 0b0u},
                      {0b1111'1110u, 0b1111'1111u, 0b1u}, 10u);
    matrix.add_column({0b0000'0001u, 0b0000'0000u, 0b0u},
                      {0b0000'1111u, 0b0000'1100u, 0b0u},
                      {0b1111'1001u, 0b1111'1111u, 0b1u}, 11u);
    matrix.add_column({0b0000'0110u, 0b0000'0000u, 0b0u},
                      {0b0011'1111u, 0b0011'0000u, 0b0u},
                      {0b1110'0011u, 0b1111'1111u, 0b1u}, 12u);
    matrix.add_column({0b0001'1100u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1100'0000u, 0b1u},
                      {0b1000'0111u, 0b1111'1111u, 0b1u}, 13u);
    matrix.add_column({0b0111'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'0011u, 0b0000'0011u, 0b0u},
                      {0b0001'1110u, 0b1111'1110u, 0b1u}, 14u);
    matrix.add_column({0b1100'0001u, 0b0000'0001u, 0b0u},
                      {0b1100'1111u, 0b0000'1111u, 0b0u},
                      {0b0111'1101u, 0b1111'1000u, 0b1u}, 15u);
    matrix.add_column({0b0000'0010u, 0b0000'0111u, 0b0u},
                      {0b0011'1111u, 0b0011'1111u, 0b0u},
                      {0b1111'0001u, 0b1110'0001u, 0b1u}, 16u);
    matrix.add_column({0b0000'1110u, 0b0001'1100u, 0b0u},
                      {0b1111'1111u, 0b1111'1100u, 0b1u},
                      {0b1100'0011u, 0b1000'0111u, 0b1u}, 17u);
    matrix.add_column({0b0000'1000u, 0b0111'0000u, 0b0u},
                      {0b1111'0011u, 0b1111'0011u, 0b1u},
                      {0b0100'1110u, 0b0001'1111u, 0b0u}, 18u);

    // row-wise matrix
    std::vector<std::vector<detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<detail::trace_directions>> expect
    {
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,NON},
        {U  ,D  ,DUL,DU ,DU ,D  ,DUL,DU ,DU ,D  },
        {U  ,DU ,D  ,DUL,DU ,DU ,D  ,DL ,DUL,DU },
        {U  ,U  ,D  ,DL ,DUL,U  ,DU ,D  ,DL ,U  },
        {U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,DL ,UL },
        {U  ,U  ,U  ,D  ,DL ,DUL,U  ,DU ,D  ,D  },
        {U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,D  },
        {U  ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,DU ,DU },
        {U  ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  },
        {NON,DU ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,DU },
        {NON,NON,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU },
        {NON,NON,NON,U  ,U  ,U  ,D  ,DL ,DUL,U  },
        {NON,NON,NON,NON,U  ,U  ,DU ,D  ,DL ,U  },
        {NON,NON,NON,NON,NON,U  ,U  ,D  ,DL ,DUL},
        {NON,NON,NON,NON,NON,NON,U  ,DU ,D  ,DL },
        {NON,NON,NON,NON,NON,NON,NON,U  ,D  ,DL },
        {NON,NON,NON,NON,NON,NON,NON,NON,DU ,D  },
        {NON,NON,NON,NON,NON,NON,NON,NON,NON,D  }
    };

    EXPECT_EQ(result, expect);
}
