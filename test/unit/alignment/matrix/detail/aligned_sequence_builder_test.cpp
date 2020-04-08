// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <list>
#include <vector>

#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/algorithm>

using seqan3::operator|;

template <typename test_type>
struct aligned_sequence_builder_fixture : ::testing::Test
{
    static constexpr seqan3::detail::trace_directions N = seqan3::detail::trace_directions::none;
    static constexpr seqan3::detail::trace_directions D = seqan3::detail::trace_directions::diagonal;
    static constexpr seqan3::detail::trace_directions U = seqan3::detail::trace_directions::up;
    static constexpr seqan3::detail::trace_directions UO = seqan3::detail::trace_directions::up_open;
    static constexpr seqan3::detail::trace_directions L = seqan3::detail::trace_directions::left;
    static constexpr seqan3::detail::trace_directions LO = seqan3::detail::trace_directions::left_open;

    seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions> matrix{seqan3::detail::number_rows{3},
                                                                                    seqan3::detail::number_cols{4},
                                                                                    std::vector
    {
        N,           LO, L,          L,
        UO, D | LO | UO, L, D | L | UO,
        U,       LO | U, D,          L
    }};

    using fst_seq_t = std::tuple_element_t<0, test_type>;
    using sec_seq_t = std::tuple_element_t<1, test_type>;

    using fst_seq_value_t = std::ranges::range_value_t<fst_seq_t>;
    using sec_seq_value_t = std::ranges::range_value_t<sec_seq_t>;

    using type_param = seqan3::detail::aligned_sequence_builder<fst_seq_t, sec_seq_t>;

    void SetUp()
    {
        fst.push_back(seqan3::assign_char_to('A', fst_seq_value_t{}));
        fst.push_back(seqan3::assign_char_to('C', fst_seq_value_t{}));
        fst.push_back(seqan3::assign_char_to('G', fst_seq_value_t{}));
        sec.push_back(seqan3::assign_char_to('A', sec_seq_value_t{}));
        sec.push_back(seqan3::assign_char_to('G', sec_seq_value_t{}));

        builder = type_param{fst, sec};
    }

    auto path(seqan3::detail::matrix_offset const & offset)
    {
        using iterator_t = decltype(seqan3::detail::trace_iterator{matrix.begin() + offset});
        return std::ranges::subrange<iterator_t, std::ranges::default_sentinel_t>
        {
            seqan3::detail::trace_iterator{matrix.begin() + offset},
            std::ranges::default_sentinel
        };
    }

    std::remove_reference_t<fst_seq_t> fst;
    std::remove_reference_t<sec_seq_t> sec;
    type_param builder;
};

using test_types = ::testing::Types<std::pair<seqan3::dna4_vector &, seqan3::dna4_vector &>,
                                    std::pair<seqan3::dna4_vector &, seqan3::dna15_vector &>,
                                    std::pair<seqan3::dna4_vector &, std::list<seqan3::dna4> &>,
                                    std::pair<std::list<seqan3::dna4> &, std::list<seqan3::dna4> &>>;

TYPED_TEST_SUITE(aligned_sequence_builder_fixture, test_types, );

TYPED_TEST(aligned_sequence_builder_fixture, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<typename TestFixture::type_param>);
    EXPECT_TRUE(std::is_copy_constructible_v<typename TestFixture::type_param>);
    EXPECT_TRUE(std::is_move_constructible_v<typename TestFixture::type_param>);
    EXPECT_TRUE(std::is_copy_assignable_v<typename TestFixture::type_param>);
    EXPECT_TRUE(std::is_move_assignable_v<typename TestFixture::type_param>);
    EXPECT_TRUE(std::is_destructible_v<typename TestFixture::type_param>);
    EXPECT_TRUE((std::is_constructible_v<typename TestFixture::type_param,
                                         typename TestFixture::fst_seq_t,
                                         typename TestFixture::sec_seq_t>));
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_2_3)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                      seqan3::detail::column_index_type{3}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0, 3}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0, 2}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"--ACG"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"AG---"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_2_2)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                      seqan3::detail::column_index_type{2}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 2u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 2u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"AC"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"AG"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_2_1)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                      seqan3::detail::column_index_type{1}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 1u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 2u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"A--"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"-AG"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_2_0)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                      seqan3::detail::column_index_type{0}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 2u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"--"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"AG"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_1_3)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                      seqan3::detail::column_index_type{3}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 3u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 1u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"ACG"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"--A"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_1_2)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                      seqan3::detail::column_index_type{2}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 2u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 1u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"-AC"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"A--"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_1_1)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                      seqan3::detail::column_index_type{1}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 1u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 1u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"A"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"A"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_1_0)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                      seqan3::detail::column_index_type{0}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 1u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"-"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"A"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_0_3)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                      seqan3::detail::column_index_type{3}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 3u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"ACG"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"---"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_0_2)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                      seqan3::detail::column_index_type{2}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 2u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"AC"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"--"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_0_1)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                      seqan3::detail::column_index_type{1}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 1u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"A"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"-"});
}

TYPED_TEST(aligned_sequence_builder_fixture, build_from_0_0)
{
    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                      seqan3::detail::column_index_type{0}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = this->builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{""});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{""});
}

TYPED_TEST(aligned_sequence_builder_fixture, both_empty)
{
    std::remove_reference_t<typename TestFixture::fst_seq_t> first{};
    std::remove_reference_t<typename TestFixture::sec_seq_t> second{};

    typename TestFixture::type_param builder{first, second};

    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                      seqan3::detail::column_index_type{0}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{""});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{""});
}

TYPED_TEST(aligned_sequence_builder_fixture, first_empty)
{
    std::remove_reference_t<typename TestFixture::fst_seq_t> first{};

    typename TestFixture::type_param builder{first, this->sec};

    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                      seqan3::detail::column_index_type{0}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 2u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"--"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"AG"});
}

TYPED_TEST(aligned_sequence_builder_fixture, second_empty)
{
    std::remove_reference_t<typename TestFixture::sec_seq_t> second{};

    typename TestFixture::type_param builder{this->fst, second};

    auto p = this->path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                      seqan3::detail::column_index_type{3}});
    auto [first_sequence_slice_positions, second_sequence_slice_positions, alignment] = builder(p);

    EXPECT_EQ(first_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 3u}));
    EXPECT_EQ(second_sequence_slice_positions, (std::pair<size_t, size_t>{0u, 0u}));
    EXPECT_EQ(std::get<0>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"ACG"});
    EXPECT_EQ(std::get<1>(alignment) | seqan3::views::to_char | seqan3::views::to<std::string>, std::string{"---"});
}
