// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/edit_distance_trace_matrix_full.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "edit_distance_trace_matrix.hpp"

TEST(semi_global_max_errors, empty)
{
    matrix_type<true, true> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, epsilon)
{
    matrix_type<true, true> matrix{1u};

    matrix.add_column({}, {}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N}};

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
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, N, N, N, N}};

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
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect
    {
        {N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  },
        {u  ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,D  },
        {u  ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,u  },
        {u  ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,u  },
        {u  ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,u  },
        {u  ,Du ,D  ,u  ,Du ,Du ,Du ,Du ,D  ,D  },
        {N  ,u  ,u  ,D  ,D  ,Dul,Du ,Du ,Du ,Du },
        {N  ,N  ,u  ,u  ,Du ,D  ,Dl ,Dul,Du ,Du },
        {N  ,N  ,N  ,u  ,Du ,u  ,D  ,D  ,Dl ,N  }
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
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect
    {
        {N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  },
        {u  ,D  ,Dul,Du ,Du ,D  ,Dul,Du ,Du ,D  },
        {u  ,Du ,D  ,Dul,Du ,Du ,D  ,Dl ,Dul,Du },
        {u  ,u  ,D  ,Dl ,Dul,u  ,Du ,D  ,Dl ,u  },
        {u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,Dl ,ul },
        {u  ,u  ,u  ,D  ,Dl ,Dul,u  ,Du ,D  ,D  },
        {u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,D  },
        {u  ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,Du ,Du },
        {u  ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  },
        {N  ,Du ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,Du },
        {N  ,N  ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du },
        {N  ,N  ,N  ,u  ,u  ,u  ,D  ,Dl ,Dul,u  },
        {N  ,N  ,N  ,N  ,u  ,u  ,Du ,D  ,Dl ,u  },
        {N  ,N  ,N  ,N  ,N  ,u  ,u  ,D  ,Dl ,Dul},
        {N  ,N  ,N  ,N  ,N  ,N  ,u  ,Du ,D  ,Dl },
        {N  ,N  ,N  ,N  ,N  ,N  ,N  ,u  ,D  ,Dl },
        {N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,Du ,D  },
        {N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,D  }
    };

    EXPECT_EQ(result, expect);
}
