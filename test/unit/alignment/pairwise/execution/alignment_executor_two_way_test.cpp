// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alignment/pairwise/execution/alignment_executor_two_way.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/view_all.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

struct dummy_alignment
{

    template <typename indexed_sequence_pairs_t>
    constexpr auto operator()(indexed_sequence_pairs_t && indexed_sequence_pairs) const
    {
        std::vector<size_t> result{};
        for (auto && indexed_pair : indexed_sequence_pairs)
        {
            using std::get;

            size_t count = 0;
            auto & [first_seq, second_seq] = get<0>(indexed_pair);
            std::ranges::for_each(views::zip(first_seq, second_seq), [&](auto && tpl)
            {
                auto && [v1, v2] = tpl;
                if (v1 == v2)
                    ++count;
            });
            result.push_back(count);
        }

        return result;
    }
};

template <typename resource_t>
struct algorithm_type_for_input
{
    using indexed_sequence_pairs_t = typename detail::chunked_indexed_sequence_pairs<resource_t>::type;
    using algorithm_input_t = std::ranges::range_value_t<indexed_sequence_pairs_t>;
    using type = std::function<std::vector<size_t>(algorithm_input_t)>;
};

template <typename t>
struct alignment_executor_two_way_test : public ::testing::Test
{
    // Some globally defined test types
    using sequence_pair_t = std::pair<std::string, std::string>;
    using sequence_pairs_t = std::vector<sequence_pair_t>;

    sequence_pair_t sequence_pair{"AACGTACGT", "ATCGTCCGT"};
    sequence_pairs_t sequence_pairs{5, sequence_pair};
};

using testing_types = testing::Types<detail::execution_handler_sequential, detail::execution_handler_parallel>;
TYPED_TEST_CASE(alignment_executor_two_way_test, testing_types);

TYPED_TEST(alignment_executor_two_way_test, construction)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<typename TestFixture::sequence_pairs_t &,
                                                                    algorithm_t,
                                                                    TypeParam>;

    EXPECT_FALSE(std::is_default_constructible_v<alignment_executor_t>);
    EXPECT_FALSE(std::is_copy_constructible_v<alignment_executor_t>);
    EXPECT_TRUE(std::is_move_constructible_v<alignment_executor_t>);
    EXPECT_FALSE(std::is_copy_assignable_v<alignment_executor_t>);
    EXPECT_TRUE(std::is_move_assignable_v<alignment_executor_t>);
}

TYPED_TEST(alignment_executor_two_way_test, construct_with_chunk_size)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<typename TestFixture::sequence_pairs_t &,
                                                                    algorithm_t,
                                                                    TypeParam>;

    {  // default chunk size should be 1.
        alignment_executor_t exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}};
        EXPECT_EQ(exec.chunk_size(), 1u);
    }

    { // chunk size should be 4.
        alignment_executor_t exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}, 4};
        EXPECT_EQ(exec.chunk_size(), 4u);
    }

    { // with chunk size and execution policy.
        alignment_executor_t exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}, 4, seq};
        EXPECT_EQ(exec.chunk_size(), 4u);
    }

    { // throw with invalid chunk size.
        EXPECT_THROW((alignment_executor_t{this->sequence_pairs, algorithm_t{dummy_alignment{}}, 0, seq}),
                     std::invalid_argument);
    }
}

TYPED_TEST(alignment_executor_two_way_test, is_eof)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<typename TestFixture::sequence_pairs_t &,
                                                                    algorithm_t,
                                                                    TypeParam>;

    alignment_executor_t exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}};
    EXPECT_FALSE(exec.is_eof());
}

TYPED_TEST(alignment_executor_two_way_test, type_deduction)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    detail::alignment_executor_two_way exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}};
    EXPECT_FALSE(exec.is_eof());
}

TYPED_TEST(alignment_executor_two_way_test, bump)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<typename TestFixture::sequence_pairs_t &,
                                                                    algorithm_t,
                                                                    TypeParam>;

    alignment_executor_t exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}};

    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, in_avail)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<typename TestFixture::sequence_pairs_t &,
                                                                    algorithm_t,
                                                                    TypeParam>;

    alignment_executor_t exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}};

    EXPECT_EQ(exec.in_avail(), 0u);
    EXPECT_EQ(exec.bump().value(), 7u);
    if constexpr (std::same_as<TypeParam, detail::execution_handler_parallel>)
        EXPECT_EQ(exec.in_avail(), 4u);
    else
        EXPECT_EQ(exec.in_avail(), 0u);
}

TYPED_TEST(alignment_executor_two_way_test, lvalue_sequence_pair_view)
{
    auto v = std::views::single(this->sequence_pair);
    using algorithm_t = typename algorithm_type_for_input<decltype(v)>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<decltype(v),
                                                                    algorithm_t,
                                                                    TypeParam>;

    alignment_executor_t exec{v, algorithm_t{dummy_alignment{}}};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, rvalue_sequence_pair_view)
{
    using single_pair_t = decltype(std::views::single(this->sequence_pair));
    using algorithm_t = typename algorithm_type_for_input<single_pair_t>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<single_pair_t,
                                                                    algorithm_t,
                                                                    TypeParam>;

    alignment_executor_t exec{std::views::single(this->sequence_pair), algorithm_t{dummy_alignment{}}};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, lvalue_sequence_pairs)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<typename TestFixture::sequence_pairs_t &,
                                                                    algorithm_t,
                                                                    TypeParam>;

    alignment_executor_t exec{this->sequence_pairs, algorithm_t{dummy_alignment{}}};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, rvalue_sequence_pairs_view)
{
    using persist_pairs_t = decltype(this->sequence_pairs | views::persist);
    using algorithm_t = typename algorithm_type_for_input<persist_pairs_t>::type;
    using alignment_executor_t = detail::alignment_executor_two_way<persist_pairs_t,
                                                                    algorithm_t,
                                                                    TypeParam>;

    alignment_executor_t exec{this->sequence_pairs | views::persist, algorithm_t{dummy_alignment{}}};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}
