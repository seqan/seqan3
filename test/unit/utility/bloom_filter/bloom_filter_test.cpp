// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/test/cereal.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

template <typename bf_type>
struct bloom_filter_test : public ::testing::Test
{
    static bf_type make_bf(seqan3::bin_size bits)
    {
        return bf_type{seqan3::bloom_filter{bits}};
    }

    static bf_type make_bf(seqan3::bin_size bits, seqan3::hash_function_count funs)
    {
        return bf_type{seqan3::bloom_filter{bits, funs}};
    }
};

using bf_types = ::testing::Types<seqan3::bloom_filter<seqan3::data_layout::uncompressed>,
                                  seqan3::bloom_filter<seqan3::data_layout::compressed>>;

TYPED_TEST_SUITE(bloom_filter_test, bf_types, );

TYPED_TEST(bloom_filter_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_copy_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_move_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_copy_assignable_v<TypeParam>);
    EXPECT_TRUE(std::is_move_assignable_v<TypeParam>);
    EXPECT_TRUE(std::is_destructible_v<TypeParam>);

    // num hash functions defaults to two
    TypeParam bf1{TestFixture::make_bf(seqan3::bin_size{1024u})};
    TypeParam bf2{TestFixture::make_bf(seqan3::bin_size{1024u}, seqan3::hash_function_count{2u})};
    EXPECT_TRUE(bf1 == bf2);

    // bin_size parameter is too small
    EXPECT_THROW((TestFixture::make_bf(seqan3::bin_size{0u})), std::logic_error);
    // not enough hash functions
    EXPECT_THROW((TestFixture::make_bf(seqan3::bin_size{32u}, seqan3::hash_function_count{0u})), std::logic_error);
    // too many hash functions
    EXPECT_THROW((TestFixture::make_bf(seqan3::bin_size{32u}, seqan3::hash_function_count{6u})), std::logic_error);
}

TYPED_TEST(bloom_filter_test, member_getter)
{
    TypeParam t1{TestFixture::make_bf(seqan3::bin_size{1024u})};
    EXPECT_EQ(t1.bit_size(), 1024u);
    EXPECT_EQ(t1.hash_function_count(), 2u);

    TypeParam t2{TestFixture::make_bf(seqan3::bin_size{1019u}, seqan3::hash_function_count{3u})};
    EXPECT_EQ(t2.bit_size(), 1019u);
    EXPECT_EQ(t2.hash_function_count(), 3u);
}

TYPED_TEST(bloom_filter_test, contains)
{
    TypeParam bf{TestFixture::make_bf(seqan3::bin_size{1024u})};

    for (size_t hash : std::views::iota(0u, 64u)) // test some hashes
    {
        // Expect false for all queries since we did not insert anything
        EXPECT_FALSE(bf.contains(hash));
    }
}

TYPED_TEST(bloom_filter_test, emplace)
{
    // 1. Test uncompressed Bloom Filter directly because the compressed one is not mutable.
    seqan3::bloom_filter bf{seqan3::bin_size{1024u}, seqan3::hash_function_count{2u}};

    for (size_t hash : std::views::iota(0u, 64u))
        bf.emplace(hash);

    // 2. Construct either the uncompressed or compressed Bloom Filter and test via contains
    TypeParam bf2{bf};
    for (size_t hash : std::views::iota(0u, 64u))
        EXPECT_TRUE(bf.contains(hash));
}

TYPED_TEST(bloom_filter_test, counting)
{
    // 1. Test uncompressed Bloom Filter directly because the compressed one is not mutable.
    seqan3::bloom_filter bf{seqan3::bin_size{1024u}, seqan3::hash_function_count{2u}};

    for (size_t hash : std::views::iota(0u, 128u))
        bf.emplace(hash);

    // 2. Construct either the uncompressed or compressed Bloom Filter and test set with bulk_contains
    TypeParam bf2{bf};

    // Test counting with all elements
    EXPECT_EQ(bf2.count(std::views::iota(0u, 128u)), 128u);

    // Test counting with some elements
    EXPECT_EQ(bf2.count(std::views::iota(22u, 42u)), 20u);
}

TYPED_TEST(bloom_filter_test, reset)
{
    // 1. Test uncompressed Bloom Filter directly because the compressed one is not mutable.
    seqan3::bloom_filter bf{seqan3::bin_size{1024u}, seqan3::hash_function_count{2u}};

    for (size_t hash : std::views::iota(0u, 64u))
        bf.emplace(hash);

    // 2. Reset the Bloom Filter
    bf.reset();

    // 3. Construct either the uncompressed or compressed Bloom Filter and test set with contains
    TypeParam bf2{bf};
    EXPECT_EQ(bf.count(std::views::iota(0u, 64u)), 0u); // nothing should be present in the Bloom Filter
}

TYPED_TEST(bloom_filter_test, data_access)
{
    seqan3::bloom_filter bf{seqan3::bin_size{1024u}};
    EXPECT_LE(sdsl::size_in_mega_bytes(bf.raw_data()), 0.001f);
}

TYPED_TEST(bloom_filter_test, serialisation)
{
    TypeParam bf{TestFixture::make_bf(seqan3::bin_size{1024u})};
    seqan3::test::do_serialisation(bf);
}
