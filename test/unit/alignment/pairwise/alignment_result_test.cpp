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

#include <type_traits>
#include <utility>

#include <meta/meta.hpp>

#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;

template <typename T>
class align_result_test : public ::testing::Test
{
public:

    using res_t   = align_result<T>;
    using id_t    = meta::at_c<T, 0>;
    using score_t = meta::at_c<T, 1>;
    using end_t   = meta::at_c<T, 2>;
    using begin_t = meta::at_c<T, 3>;
    using trace_t = meta::at_c<T, 4>;
    using seq_t   = std::tuple_element_t<0, trace_t>;
};

using aligned_seq = std::vector<gapped<dna4>>;
// add all alphabets here
using align_result_test_types = ::testing::Types<type_list<uint32_t,
                                                 int32_t,
                                                 std::pair<size_t, size_t>,
                                                 std::pair<size_t, size_t>,
                                                 std::tuple<aligned_seq, aligned_seq>>>;

TYPED_TEST_CASE(align_result_test, align_result_test_types);

TYPED_TEST(align_result_test, constructor)
{
    using res_t = typename TestFixture::res_t;
    EXPECT_TRUE((std::is_default_constructible_v<res_t>));
    EXPECT_TRUE((std::is_copy_constructible_v<res_t>));
    EXPECT_TRUE((std::is_move_constructible_v<res_t>));
    EXPECT_TRUE((std::is_destructible_v<res_t>));
}

TYPED_TEST(align_result_test, tuple_concept)
{
    EXPECT_TRUE((tuple_like_concept<typename TestFixture::res_t>));
}

TYPED_TEST(align_result_test, tuple_element)
{
    using res_t = typename TestFixture::res_t;
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, res_t>, typename TestFixture::id_t>));
    EXPECT_TRUE((std::is_same_v<std::tuple_element_t<4, res_t>, typename TestFixture::trace_t>));
}

TYPED_TEST(align_result_test, tuple_size)
{
    using res_t = typename TestFixture::res_t;
    EXPECT_EQ(std::tuple_size_v<res_t>, 5u);
}

TYPED_TEST(align_result_test, std_position_get)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(std::get<0>(tmp), 1u);
        EXPECT_EQ(std::string{std::get<0>(std::get<4>(tmp)) | view::to_char}, std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(std::get<0>(tmp)), typename TestFixture::id_t &>));
    }

    {
        res_t const tmp{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(std::get<0>(tmp), 1u);
        EXPECT_EQ(std::string{std::get<0>(std::get<4>(tmp)) | view::to_char}, std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(std::get<0>(tmp)), typename TestFixture::id_t const &>));
    }

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(std::get<0>(std::move(tmp)), 1u);
        res_t tmp2{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(std::string{view::persist(std::get<0>(std::get<4>(std::move(tmp2)))) | view::to_char},
                  std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(std::get<0>(std::move(tmp))), typename TestFixture::id_t &&>));
    }

    {
        res_t const tmp{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(std::get<0>(std::move(tmp)), 1u);
        res_t const tmp2{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(std::string{view::persist(std::get<0>(std::get<4>(std::move(tmp2)))) | view::to_char},
                  std::string{"AT-C--A"});
        // TODO: enable if gcc-7 is fixed.
        // EXPECT_TRUE((std::is_same_v<decltype(std::get<0>(std::move(tmp))), typename TestFixture::id_t const &&>));
    }
}

TYPED_TEST(align_result_test, seqan3_pos_get)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(get<0>(tmp), 1u);
        EXPECT_EQ(std::string{std::get<0>(get<4>(tmp)) | view::to_char}, std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(get<0>(tmp)), typename TestFixture::id_t &>));
    }

    {
        res_t const tmp{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(get<0>(tmp), 1u);
        EXPECT_EQ(std::string{std::get<0>(get<4>(tmp)) | view::to_char}, std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(get<0>(tmp)), typename TestFixture::id_t const &>));
    }

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(get<0>(std::move(tmp)), 1u);
        res_t tmp2{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(std::string{view::persist(std::get<0>(get<4>(std::move(tmp2)))) | view::to_char},
                  std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(get<0>(std::move(tmp))), typename TestFixture::id_t &&>));
    }

    {
        res_t const tmp{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(get<0>(std::move(tmp)), 1u);
        res_t const tmp2{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(std::string{view::persist(std::get<0>(get<4>(std::move(tmp2)))) | view::to_char},
                  std::string{"AT-C--A"});
        // TODO: enable if gcc-7 is fixed.
        // EXPECT_TRUE((std::is_same_v<decltype(get<0>(std::move(tmp))), typename TestFixture::id_t const &&>));
    }
}

TYPED_TEST(align_result_test, seqan3_enum_get)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(get<align_result_key::id>(tmp), 1u);
        EXPECT_EQ(std::string{std::get<0>(get<align_result_key::trace>(tmp)) | view::to_char}, std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(get<align_result_key::id>(tmp)), typename TestFixture::id_t &>));
    }

    {
        res_t const tmp{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(get<align_result_key::id>(tmp), 1u);
        EXPECT_EQ(std::string{std::get<0>(get<align_result_key::trace>(tmp)) | view::to_char}, std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(get<align_result_key::id>(tmp)), typename TestFixture::id_t const &>));
    }

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(get<align_result_key::id>(std::move(tmp)), 1u);
        res_t tmp2{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(std::string{view::persist(std::get<0>(get<align_result_key::trace>(std::move(tmp2)))) | view::to_char},
                  std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(get<align_result_key::id>(std::move(tmp))), typename TestFixture::id_t &&>));
    }

    {
        res_t const tmp{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(get<align_result_key::id>(std::move(tmp)), 1u);
        res_t const tmp2{res_t{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}}};
        EXPECT_EQ(std::string{view::persist(std::get<0>(get<align_result_key::trace>(std::move(tmp2)))) | view::to_char},
                  std::string{"AT-C--A"});
        EXPECT_TRUE((std::is_same_v<decltype(get<align_result_key::id>(std::move(tmp))), typename TestFixture::id_t const &&>));
    }
}

TYPED_TEST(align_result_test, id)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
    EXPECT_EQ(tmp.id(), 1u);
}

TYPED_TEST(align_result_test, score)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
    EXPECT_EQ(tmp.score(), 0);
}

TYPED_TEST(align_result_test, end_coordinate)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(tmp.end_coordinate(), (std::pair{10lu, 10lu}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.end_coordinate()),
                                    std::pair<size_t, size_t> const &>));
    }

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_TRUE(std::move(tmp).end_coordinate() == (std::pair{10lu, 10lu}));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).end_coordinate()),
                                    std::pair<size_t, size_t> const &&>));
    }

    {
        align_result<type_list<uint32_t, int32_t>> tmp{1, 0};
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(tmp.end_coordinate())>,
                                    remove_cvref_t<decltype(std::ignore)>>));
    }
}

TYPED_TEST(align_result_test, begin_coordinate)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(tmp.begin_coordinate(), (std::pair{0lu, 0lu}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.begin_coordinate()),
                                    std::pair<size_t, size_t> const &>));
    }

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_TRUE(std::move(tmp).begin_coordinate() == (std::pair{0lu, 0lu}));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).begin_coordinate()),
                                    std::pair<size_t, size_t> const &&>));
    }

    {
        align_result<type_list<uint32_t, int32_t>> tmp{1, 0};
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(tmp.begin_coordinate())>,
                                    remove_cvref_t<decltype(std::ignore)>>));
    }
}

TYPED_TEST(align_result_test, trace)
{
    using res_t = typename TestFixture::res_t;
    using seq_t = typename TestFixture::seq_t;
    seq_t seq{'A'_dna4, 'T'_dna4, gap::GAP, 'C'_dna4, gap::GAP, gap::GAP, 'A'_dna4};

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_EQ(tmp.trace(), (std::tuple{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.trace()),
                                    std::tuple<seq_t, seq_t> const &>));
    }

    {
        res_t tmp{1u, 0, {10u, 10u}, {0u, 0u}, {seq, seq}};
        EXPECT_TRUE(std::move(tmp).trace() == (std::tuple{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).trace()),
                                    std::tuple<seq_t, seq_t> const &&>));
    }

    {
        align_result<type_list<uint32_t, int32_t>> tmp{1, 0};
        EXPECT_TRUE((std::is_same_v<remove_cvref_t<decltype(tmp.trace())>,
                                    remove_cvref_t<decltype(std::ignore)>>));
    }
}
