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

#include <gtest/gtest.h>

#include <meta/meta.hpp>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;

template <typename T>
class align_result_test : public ::testing::Test
{
public:

    using id_t     = meta::at_c<T, 0>;
    using score_t  = meta::at_c<T, 1>;
    using pos_t    = meta::at_c<T, 2>;
    using trace_t  = meta::at_c<T, 3>;
    using seq_t    = meta::at_c<T, 4>;
    using coord_t  = std::pair<pos_t, pos_t>;
    using traits_t = detail::align_result_value_type<id_t, score_t, coord_t, trace_t>;
    using res_t    = align_result<traits_t>;
};

using aligned_seq = std::vector<gapped<dna4>>;
using align_result_test_types = ::testing::Types
    <type_list<uint32_t, int32_t, size_t, std::tuple<aligned_seq, aligned_seq>, aligned_seq>,
     type_list<uint32_t, int32_t, size_t, std::pair<aligned_seq, aligned_seq>,  aligned_seq>,
     type_list<uint32_t, int32_t, size_t, std::vector<aligned_seq>,             aligned_seq>,
     type_list<uint32_t, float,   size_t, std::tuple<aligned_seq, aligned_seq>, aligned_seq>,
     type_list<uint32_t, float,   size_t, std::pair<aligned_seq, aligned_seq>,  aligned_seq>,
     type_list<uint32_t, float,   size_t, std::vector<aligned_seq>,             aligned_seq>>;

TYPED_TEST_CASE(align_result_test, align_result_test_types);

TYPED_TEST(align_result_test, constructor)
{
    using res_t = typename TestFixture::res_t;
    EXPECT_TRUE((std::is_default_constructible_v<res_t>));
    EXPECT_TRUE((std::is_copy_constructible_v<res_t>));
    EXPECT_TRUE((std::is_move_constructible_v<res_t>));
    EXPECT_TRUE((std::is_copy_assignable_v<res_t>));
    EXPECT_TRUE((std::is_move_assignable_v<res_t>));
    EXPECT_TRUE((std::is_destructible_v<res_t>));
}

TYPED_TEST(align_result_test, get_id)
{
    using traits_t = typename TestFixture::traits_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        align_result<traits_t> tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_id()), typename TestFixture::id_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_id()), typename TestFixture::id_t>));
        EXPECT_EQ(tmp.get_id(), 1u);
    }

    {
        align_result<traits_t> const tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_id(), 1u);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_id()), typename TestFixture::id_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_id()), typename TestFixture::id_t>));
    }
}

TYPED_TEST(align_result_test, get_score)
{
    using traits_t = typename TestFixture::traits_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        align_result<traits_t> tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_score()), typename TestFixture::score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_score()), typename TestFixture::score_t>));
    }

    {
        align_result<traits_t> const tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_score()), typename TestFixture::score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_score()), typename TestFixture::score_t>));
    }

}

TYPED_TEST(align_result_test, end_coordinate)
{
    using traits_t = typename TestFixture::traits_t;
    using seq_t = typename TestFixture::seq_t;
    using coord_t = typename TestFixture::coord_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        align_result<traits_t> tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_end_coordinate(), (coord_t{10ul, 10ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_end_coordinate()), coord_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_end_coordinate()), coord_t>));
    }

    {
        align_result<traits_t> const tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_end_coordinate(), (coord_t{10ul, 10ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_end_coordinate()), coord_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_end_coordinate()), coord_t>));
    }
}

TYPED_TEST(align_result_test, begin_coordinate)
{
    using traits_t = typename TestFixture::traits_t;
    using seq_t = typename TestFixture::seq_t;
    using coord_t = typename TestFixture::coord_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        align_result<traits_t> tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_begin_coordinate(), (coord_t{0ul, 0ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_begin_coordinate()), coord_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_begin_coordinate()), coord_t>));
    }

    {
        align_result<traits_t> const tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_begin_coordinate(), (coord_t{0ul, 0ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_begin_coordinate()), coord_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_begin_coordinate()), coord_t>));
    }
}

TYPED_TEST(align_result_test, trace)
{
    using traits_t = typename TestFixture::traits_t;
    using seq_t = typename TestFixture::seq_t;
    using trace_t = typename TestFixture::trace_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        align_result<traits_t> tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_trace(), (trace_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_trace()), trace_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_trace()), trace_t>));
    }

    {
        align_result<traits_t> const tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.get_trace(), (trace_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_trace()), trace_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).get_trace()), trace_t>));
    }

    if constexpr (tuple_like_concept<trace_t>)
    {
        align_result<traits_t> tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(std::string{std::get<0>(tmp.get_trace()) | view::persist | view::to_char}, std::string{"AT-C--A"});
        EXPECT_EQ(std::string{std::get<1>(tmp.get_trace()) | view::persist | view::to_char}, std::string{"AT-C--A"});
    }
    else
    {
        align_result<traits_t> tmp{traits_t{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(std::string{(tmp.get_trace()[0]) | view::persist | view::to_char}, std::string{"AT-C--A"});
        EXPECT_EQ(std::string{(tmp.get_trace()[1]) | view::persist | view::to_char}, std::string{"AT-C--A"});
    }
}

TEST(align_result_test, reduced_type)
{
    {
        detail::align_result_value_type tr{2u, 5};
        align_result tmp(tr);
        EXPECT_EQ(tmp.get_id(), 2u);
        EXPECT_EQ(tmp.get_score(), 5);
    }

    {
        detail::align_result_value_type tr{2, 5.0f, std::pair<int, int>{1, -1}};
        align_result tmp(tr);
        EXPECT_EQ(tmp.get_id(), 2);
        EXPECT_FLOAT_EQ(tmp.get_score(), 5.0f);
        EXPECT_EQ((tmp.get_end_coordinate()), (std::pair<int, int>{1, -1}));
    }

    {
        detail::align_result_value_type tr{2, 5.0f, std::pair<int, int>{1, -1}, std::pair<int, int>{10, -10}};
        align_result tmp(tr);
        EXPECT_EQ(tmp.get_id(), 2);
        EXPECT_FLOAT_EQ(tmp.get_score(), 5.0f);
        EXPECT_EQ((tmp.get_end_coordinate()), (std::pair<int, int>{1, -1}));
        EXPECT_EQ((tmp.get_begin_coordinate()), (std::pair<int, int>{10, -10}));
    }
}

TEST(align_result_test, type_deduction)
{
    {
        using coord_t = std::pair<int, int>;
        std::vector<gapped<rna5>> seq{'A'_rna5, 'U'_rna5, gap{}, 'C'_rna5, gap{}, gap{}, 'A'_rna5};

        detail::align_result_value_type tr{2, 5.0, coord_t{1, -1}, coord_t{10,-10}, seq};
        align_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_end_coordinate()), coord_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_begin_coordinate()), coord_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_trace()), std::vector<gapped<rna5>>>));

        EXPECT_EQ(tmp.get_id(), 2);
        EXPECT_DOUBLE_EQ(tmp.get_score(), 5.0);
        EXPECT_EQ((tmp.get_end_coordinate()), (coord_t{1, -1}));
        EXPECT_EQ((tmp.get_begin_coordinate()), (coord_t{10, -10}));
        EXPECT_EQ(tmp.get_trace(), seq);
    }

    {
        using coord_t = std::pair<int, int>;

        detail::align_result_value_type tr{2, 5.0, coord_t{1, -1}};
        align_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_end_coordinate()), coord_t>));
    }

    {
        detail::align_result_value_type tr{2, 5.0};
        align_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.get_score()), double>));
    }
}
