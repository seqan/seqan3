// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/cereal.hpp>
#include <seqan3/test/expect_range_eq.hpp>

template <typename ibf_type>
struct interleaved_bloom_filter_test : public ::testing::Test
{
    static ibf_type make_ibf(seqan3::bin_count bins, seqan3::bin_size bits)
    {
        return ibf_type{seqan3::interleaved_bloom_filter{bins, bits}};
    }

    static ibf_type make_ibf(seqan3::bin_count bins, seqan3::bin_size bits, seqan3::hash_function_count funs)
    {
        return ibf_type{seqan3::interleaved_bloom_filter{bins, bits, funs}};
    }
};

using ibf_types = ::testing::Types<seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>,
                                   seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>>;

TYPED_TEST_SUITE(interleaved_bloom_filter_test, ibf_types, );

TYPED_TEST(interleaved_bloom_filter_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_copy_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_move_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_copy_assignable_v<TypeParam>);
    EXPECT_TRUE(std::is_move_assignable_v<TypeParam>);
    EXPECT_TRUE(std::is_destructible_v<TypeParam>);

    // num hash functions defaults to two
    TypeParam ibf1{TestFixture::make_ibf(seqan3::bin_count{64u}, seqan3::bin_size{1024u})};
    TypeParam ibf2{
        TestFixture::make_ibf(seqan3::bin_count{64u}, seqan3::bin_size{1024u}, seqan3::hash_function_count{2u})};
    EXPECT_TRUE(ibf1 == ibf2);

    // bin_size parameter is too small
    EXPECT_THROW((TestFixture::make_ibf(seqan3::bin_count{64u}, seqan3::bin_size{0u})), std::logic_error);
    // not enough bins
    EXPECT_THROW((TestFixture::make_ibf(seqan3::bin_count{0u}, seqan3::bin_size{32u})), std::logic_error);
    // not enough hash functions
    EXPECT_THROW(
        (TestFixture::make_ibf(seqan3::bin_count{64u}, seqan3::bin_size{32u}, seqan3::hash_function_count{0u})),
        std::logic_error);
    // too many hash functions
    EXPECT_THROW(
        (TestFixture::make_ibf(seqan3::bin_count{64u}, seqan3::bin_size{32u}, seqan3::hash_function_count{6u})),
        std::logic_error);
}

TYPED_TEST(interleaved_bloom_filter_test, member_getter)
{
    TypeParam t1{TestFixture::make_ibf(seqan3::bin_count{64u}, seqan3::bin_size{1024u})};
    EXPECT_EQ(t1.bin_count(), 64u);
    EXPECT_EQ(t1.bin_size(), 1024u);
    EXPECT_EQ(t1.bit_size(), 65'536ull);
    EXPECT_EQ(t1.hash_function_count(), 2u);

    TypeParam t2{
        TestFixture::make_ibf(seqan3::bin_count{73u}, seqan3::bin_size{1019u}, seqan3::hash_function_count{3u})};
    EXPECT_EQ(t2.bin_count(), 73u);
    EXPECT_EQ(t2.bin_size(), 1019u);
    EXPECT_EQ(t2.bit_size(), 130'432ull);
    EXPECT_EQ(t2.hash_function_count(), 3u);
}

TYPED_TEST(interleaved_bloom_filter_test, bulk_contains)
{
    TypeParam ibf{TestFixture::make_ibf(seqan3::bin_count{64u}, seqan3::bin_size{1024u})};
    std::vector<bool> expected(64); // empty bitvector is expected since we did not insert anything
    auto agent = ibf.membership_agent();

    for (size_t hash : std::views::iota(0, 64)) // test correct resize for each bin individually
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }

    // Test iterator interface
    for (size_t hash : std::views::iota(0, 64)) // test correct resize for each bin individually
    {
        auto & res = agent.bulk_contains(hash);
        size_t i = 0;
        for (auto it = res.begin(); it < res.end(); ++it, ++i)
        {
            EXPECT_EQ(*it, expected[i]);
        }
        EXPECT_EQ(i, expected.size());
    }

    // Test operator[] interface
    for (size_t hash : std::views::iota(0, 64)) // test correct resize for each bin individually
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_EQ(expected.size(), res.size());
        for (size_t i = 0; i < res.size(); ++i)
        {
            EXPECT_EQ(res[i], expected[i]);
        }
    }
}

TYPED_TEST(interleaved_bloom_filter_test, emplace)
{
    // 1. Test uncompressed interleaved_bloom_filter directly because the compressed one is not mutable.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{64u},
                                         seqan3::bin_size{1024u},
                                         seqan3::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    // 2. Construct either the uncompressed or compressed interleaved_bloom_filter and test set with bulk_contains
    TypeParam ibf2{ibf};
    auto agent = ibf2.membership_agent();
    std::vector<bool> expected(64, 1);          // every hash value should be set for every bin
    for (size_t hash : std::views::iota(0, 64)) // test correct resize for each bin individually
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }
}

TYPED_TEST(interleaved_bloom_filter_test, clear)
{
    // 1. Test uncompressed interleaved_bloom_filter directly because the compressed one is not mutable.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{64u},
                                         seqan3::bin_size{1024u},
                                         seqan3::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    // 2. Clear a bin
    ibf.clear(seqan3::bin_index{17u});

    // 3. Construct either the uncompressed or compressed interleaved_bloom_filter and test set with bulk_contains
    TypeParam ibf2{ibf};
    auto agent = ibf2.membership_agent();
    std::vector<bool> expected(64, 1); // every hash value should be set for every bin...
    expected[17] = 0;                  // ...except bin 17
    for (size_t hash : std::views::iota(0, 64))
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }
}

TYPED_TEST(interleaved_bloom_filter_test, clear_range)
{
    // 1. Test uncompressed interleaved_bloom_filter directly because the compressed one is not mutable.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{64u},
                                         seqan3::bin_size{1024u},
                                         seqan3::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    // 2. Clear a range of bins
    std::vector<seqan3::bin_index> bin_range{seqan3::bin_index{8u}, seqan3::bin_index{17u}, seqan3::bin_index{45u}};
    ibf.clear(bin_range);

    // 3. Construct either the uncompressed or compressed interleaved_bloom_filter and test set with bulk_contains
    TypeParam ibf2{ibf};
    auto agent = ibf2.membership_agent();
    std::vector<bool> expected(64, 1); // every hash value should be set for every bin...
    expected[8] = 0;                   // ...except bin 8
    expected[17] = 0;                  // ...except bin 17
    expected[45] = 0;                  // ...except bin 45
    for (size_t hash : std::views::iota(0, 64))
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }
}

TYPED_TEST(interleaved_bloom_filter_test, counting)
{
    // 1. Test uncompressed interleaved_bloom_filter directly because the compressed one is not mutable.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{128u},
                                         seqan3::bin_size{1024u},
                                         seqan3::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 128))
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    // 2. Construct either the uncompressed or compressed interleaved_bloom_filter and test set with bulk_contains
    TypeParam ibf2{ibf};
    seqan3::counting_vector<size_t> counting(128, 0);
    auto agent = ibf2.membership_agent();
    for (size_t hash : std::views::iota(0, 128)) // test correct resize for each bin individually
    {
        counting += agent.bulk_contains(hash);
    }
    std::vector<size_t> expected(128, 128);
    EXPECT_EQ(counting, expected);

    // Counting vectors can be added together
    std::vector<size_t> expected2(128, 256);
    counting += counting;
    EXPECT_EQ(counting, expected2);

    // minus binning_bitvector
    for (size_t hash : std::views::iota(0, 128)) // test correct resize for each bin individually
    {
        counting -= agent.bulk_contains(hash);
    }
    EXPECT_EQ(counting, std::vector<size_t>(128, 128));

    // minus other counting vector
    counting -= seqan3::counting_vector<size_t>(128, 128 - 42);
    EXPECT_EQ(counting, std::vector<size_t>(128, 42));
}

TYPED_TEST(interleaved_bloom_filter_test, counting_agent)
{
    // 1. Test uncompressed interleaved_bloom_filter directly because the compressed one is not mutable.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{128u},
                                         seqan3::bin_size{1024u},
                                         seqan3::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 128))
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    // 2. Construct either the uncompressed or compressed interleaved_bloom_filter and test set with bulk_count
    TypeParam ibf2{ibf};
    auto agent = ibf2.counting_agent();
    auto agent2 = ibf2.template counting_agent<size_t>();

    std::vector<size_t> expected(128, 128);
    EXPECT_RANGE_EQ(agent.bulk_count(std::views::iota(0u, 128u)), expected);
    EXPECT_RANGE_EQ(agent2.bulk_count(std::views::iota(0u, 128u)), expected);
}

// Check special case where there is only one `1` in the bitvector.
TYPED_TEST(interleaved_bloom_filter_test, counting_no_ub)
{
    // 1. Test uncompressed interleaved_bloom_filter directly because the compressed one is not mutable.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{128u},
                                         seqan3::bin_size{1024u},
                                         seqan3::hash_function_count{2u}};

    for (size_t bin_idx : std::array<size_t, 2>{63, 127})
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    // 2. Construct either the uncompressed or compressed interleaved_bloom_filter and test set with bulk_contains
    TypeParam ibf2{ibf};
    seqan3::counting_vector<size_t> counting(128, 0);
    auto agent = ibf2.membership_agent();
    for (size_t hash : std::views::iota(0, 128)) // test correct resize for each bin individually
    {
        counting += agent.bulk_contains(hash);
    }
    std::vector<size_t> expected(128, 0);
    expected[63] = 128;
    expected[127] = 128;
    EXPECT_EQ(counting, expected);

    // Counting vectors can be added together
    std::vector<size_t> expected2(128, 0);
    expected2[63] = 256;
    expected2[127] = 256;
    counting += counting;
    EXPECT_EQ(counting, expected2);
}

// Check special case where there is only one `1` in the bitvector.
TYPED_TEST(interleaved_bloom_filter_test, counting_agent_no_ub)
{
    // 1. Test uncompressed interleaved_bloom_filter directly because the compressed one is not mutable.
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{128u},
                                         seqan3::bin_size{1024u},
                                         seqan3::hash_function_count{2u}};

    for (size_t bin_idx : std::array<size_t, 2>{63, 127})
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    // 2. Construct either the uncompressed or compressed interleaved_bloom_filter and test set with bulk_contains
    TypeParam ibf2{ibf};
    auto agent = ibf2.counting_agent();
    auto agent2 = ibf2.template counting_agent<size_t>();

    std::vector<size_t> expected(128, 0);
    expected[63] = 128;
    expected[127] = 128;
    EXPECT_RANGE_EQ(agent.bulk_count(std::views::iota(0u, 128u)), expected);
    EXPECT_RANGE_EQ(agent2.bulk_count(std::views::iota(0u, 128u)), expected);
}

TYPED_TEST(interleaved_bloom_filter_test, increase_bin_number_to)
{
    seqan3::interleaved_bloom_filter ibf1{seqan3::bin_count{73u}, seqan3::bin_size{1024u}};
    seqan3::interleaved_bloom_filter ibf2{ibf1};

    // 1. Throw if trying to reduce number of bins.
    EXPECT_THROW(ibf1.increase_bin_number_to(seqan3::bin_count{62u}), std::invalid_argument);

    // 2. No change in bin_words implies no change in size.
    ibf2.increase_bin_number_to({seqan3::bin_count{127u}});
    EXPECT_EQ(ibf1.bit_size(), ibf2.bit_size());
    EXPECT_EQ(ibf2.bin_count(), 127u);

    // 3. If resizing takes place, the inserted values must still be valid.
    auto hashes = std::views::iota(0, 64);
    for (size_t current_bin : std::views::iota(0, 64)) // test correct resize for each bin individually
    {
        seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{64u}, seqan3::bin_size{1024u}};
        std::ranges::for_each(hashes,
                              [&ibf, &current_bin](auto const h)
                              {
                                  ibf.emplace(h, seqan3::bin_index{current_bin});
                              });

        ibf.increase_bin_number_to(seqan3::bin_count{73u});

        EXPECT_EQ(ibf.bin_count(), 73u);
        EXPECT_GE(ibf.bit_size(), 1024u);

        std::vector<bool> expected(73, 0);
        expected[current_bin] = 1; // none of the bins except current_bin stores the hash values.
        TypeParam tibf{ibf};       // test output on compressed and uncompressed
        auto agent = tibf.membership_agent();
        for (size_t const h : hashes)
        {
            auto & res = agent.bulk_contains(h);
            EXPECT_RANGE_EQ(res, expected);
        }
    }
}

TYPED_TEST(interleaved_bloom_filter_test, data_access)
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{1024u}, seqan3::bin_size{1024u}};

    EXPECT_LE(sdsl::size_in_mega_bytes(ibf.raw_data()), 1.0f);
}

TYPED_TEST(interleaved_bloom_filter_test, serialisation)
{
    TypeParam ibf{TestFixture::make_ibf(seqan3::bin_count{73u}, seqan3::bin_size{1024u})};
    seqan3::test::do_serialisation(ibf);
}

TEST(interleaved_bloom_filter_test, decompression)
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{64u}, seqan3::bin_size{1024u}};

    // Only use every other bin.
    auto take_odd = [](auto number)
    {
        return number & 1;
    };
    for (size_t bin_idx : std::views::iota(0, 64) | std::views::filter(take_odd))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan3::bin_index{bin_idx});

    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_compressed{ibf};

    seqan3::interleaved_bloom_filter ibf_decompressed{ibf_compressed};

    EXPECT_TRUE(ibf == ibf_decompressed);
}
