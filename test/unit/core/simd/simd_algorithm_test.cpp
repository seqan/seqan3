// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/std/span>
#include <seqan3/test/simd_utility.hpp>

using namespace seqan3;

TEST(simd_algorithm, fill)
{
    using simd_type = simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < simd_traits<simd_type>::length; ++i)
        expect[i] = 4;

    constexpr simd_type result = fill<simd_type>(4);
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, iota)
{
    using simd_type = simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < simd_traits<simd_type>::length; ++i)
        expect[i] = i;

    constexpr simd_type result = iota<simd_type>(0);
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, transform_to_soa)
{
    using simd_t = simd_type_t<int32_t, 4>;

    EXPECT_EQ(simd_traits<simd_t>::length, 4);

    auto seq1 = "ATGCAAAAATA"_dna4;
    auto seq2 = "CATGCCCCCGC"_dna4;
    auto seq3 = "GCATGGGGGGC"_dna4;
    auto seq4 = "TGCATTTTTTA"_dna4;

    std::vector data{std::span{std::ranges::data(seq1), simd_traits<simd_t>::length},
                     std::span{std::ranges::data(seq2), simd_traits<simd_t>::length},
                     std::span{std::ranges::data(seq3), simd_traits<simd_t>::length},
                     std::span{std::ranges::data(seq4), simd_traits<simd_t>::length}};

    std::vector<simd_t, aligned_allocator<simd_t, simd_traits<simd_t>::max_length>> out_vec;

    detail::transform_batch_to_soa<simd_t>(std::back_inserter(out_vec), data);

    SIMD_EQ(out_vec[0], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[1], simd_t{3, 0, 1, 2});
    SIMD_EQ(out_vec[2], simd_t{2, 3, 0, 1});
    SIMD_EQ(out_vec[3], simd_t{1, 2, 3, 0});

    data[0] = std::span{std::ranges::data(seq1) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};
    data[1] = std::span{std::ranges::data(seq2) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};
    data[2] = std::span{std::ranges::data(seq3) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};
    data[3] = std::span{std::ranges::data(seq4) + simd_traits<simd_t>::length, simd_traits<simd_t>::length};

    detail::transform_batch_to_soa<simd_t>(std::back_inserter(out_vec), data);

    SIMD_EQ(out_vec[4], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[5], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[6], simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[7], simd_t{0, 1, 2, 3});

    data[0] = std::span{std::ranges::data(seq1) + 2 * simd_traits<simd_t>::length, 2};
    data[1] = std::span{std::ranges::data(seq2) + 2 * simd_traits<simd_t>::length, 2};
    data[2] = std::span{std::ranges::data(seq3) + 2 * simd_traits<simd_t>::length, 2};
    data[3] = std::span{std::ranges::data(seq4) + 2 * simd_traits<simd_t>::length, 2};

    detail::transform_batch_to_soa<simd_t>(std::back_inserter(out_vec), data);

    SIMD_EQ(out_vec[8],  simd_t{0, 1, 2, 3});
    SIMD_EQ(out_vec[9],  simd_t{3, 2, 2, 3});
    SIMD_EQ(out_vec[10], simd_t{0, 1, 1, 0});
}
