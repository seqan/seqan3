// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <numeric>

#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/test/simd_utility.hpp>

TEST(simd_algorithm, fill)
{
    using simd_type = seqan3::simd::simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < seqan3::simd::simd_traits<simd_type>::length; ++i)
        expect[i] = 4;

    constexpr simd_type result = seqan3::simd::fill<simd_type>(4);
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, iota)
{
    using simd_type = seqan3::simd::simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < seqan3::simd::simd_traits<simd_type>::length; ++i)
        expect[i] = i;

    constexpr simd_type result = seqan3::simd::iota<simd_type>(0);
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, transpose)
{
    using simd_t = seqan3::simd::simd_type_t<uint8_t>;

    if constexpr (seqan3::simd::simd_traits<simd_t>::length > 1)
    {
        std::array<simd_t, seqan3::simd::simd_traits<simd_t>::length> matrix;

        for (size_t i = 0; i < matrix.size(); ++i)
            matrix[i] = seqan3::simd::iota<simd_t>(0);

        seqan3::simd::transpose(matrix);

        for (size_t i = 0; i < matrix.size(); ++i)
            SIMD_EQ(matrix[i], seqan3::simd::fill<simd_t>(i));
    }
}

//-----------------------------------------------------------------------------
// Algorithm load
//-----------------------------------------------------------------------------

template <typename simd_t>
struct simd_algorithm_load : ::testing::Test
{

    void SetUp()
    {
        memory.resize(100);
        std::iota(memory.begin(), memory.end(), 0);
    }

    std::vector<typename seqan3::simd::simd_traits<simd_t>::scalar_type> memory;
};

using simd_load_types = ::testing::Types<seqan3::simd::simd_type_t<int8_t>,
                                         seqan3::simd::simd_type_t<uint8_t>,
                                         seqan3::simd::simd_type_t<int16_t>,
                                         seqan3::simd::simd_type_t<uint16_t>,
                                         seqan3::simd::simd_type_t<int32_t>,
                                         seqan3::simd::simd_type_t<uint32_t>,
                                         seqan3::simd::simd_type_t<int64_t>,
                                         seqan3::simd::simd_type_t<uint64_t>>;

TYPED_TEST_SUITE(simd_algorithm_load, simd_load_types, );

TYPED_TEST(simd_algorithm_load, load)
{
    SIMD_EQ(seqan3::simd::load<TypeParam>(this->memory.data()), seqan3::simd::iota<TypeParam>(0));
    SIMD_EQ(seqan3::simd::load<TypeParam>(this->memory.data() + 10), seqan3::simd::iota<TypeParam>(10));
}

//-----------------------------------------------------------------------------
// Algorithm extract
//-----------------------------------------------------------------------------

template <typename simd_t>
struct simd_algorithm_extract : ::testing::Test
{
    static constexpr size_t simd_length = seqan3::simd::simd_traits<simd_t>::length;
};

using simd_extract_types = ::testing::Types<seqan3::simd::simd_type_t<uint8_t>,
                                            seqan3::simd::simd_type_t<uint16_t>,
                                            seqan3::simd::simd_type_t<int32_t>,
                                            seqan3::simd::simd_type_t<int64_t>>;
TYPED_TEST_SUITE(simd_algorithm_extract, simd_extract_types, );

TYPED_TEST(simd_algorithm_extract, extract_halve)
{
    TypeParam vec = seqan3::simd::iota<TypeParam>(0);

    // + 1 needed for emulated types without arch specification (simd length = 1).
    for (size_t idx = 0; idx < (TestFixture::simd_length + 1)/ 2; ++idx)
    {
        EXPECT_EQ(seqan3::detail::extract_halve<0>(vec)[idx], vec[idx]);
        EXPECT_EQ(seqan3::detail::extract_halve<1>(vec)[idx], vec[idx + TestFixture::simd_length / 2]);
    }
}

TYPED_TEST(simd_algorithm_extract, extract_quarter)
{
    TypeParam vec = seqan3::simd::iota<TypeParam>(0);

    // + 1 needed for emulated types without arch specification (simd length = 1).
    for (size_t idx = 0; idx < (TestFixture::simd_length + 1) / 4; ++idx)
    {
        EXPECT_EQ(seqan3::detail::extract_quarter<0>(vec)[idx], vec[idx + TestFixture::simd_length / 4 * 0]);
        EXPECT_EQ(seqan3::detail::extract_quarter<1>(vec)[idx], vec[idx + TestFixture::simd_length / 4 * 1]);
        EXPECT_EQ(seqan3::detail::extract_quarter<2>(vec)[idx], vec[idx + TestFixture::simd_length / 4 * 2]);
        EXPECT_EQ(seqan3::detail::extract_quarter<3>(vec)[idx], vec[idx + TestFixture::simd_length / 4 * 3]);
    }
}

TYPED_TEST(simd_algorithm_extract, extract_eighth)
{
    TypeParam vec = seqan3::simd::iota<TypeParam>(0);

    // + 1 needed for emulated types without arch specification (simd length = 1).
    for (size_t idx = 0; idx < (TestFixture::simd_length + 1) / 8; ++idx)
    {
        EXPECT_EQ(seqan3::detail::extract_eighth<0>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 0]);
        EXPECT_EQ(seqan3::detail::extract_eighth<1>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 1]);
        EXPECT_EQ(seqan3::detail::extract_eighth<2>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 2]);
        EXPECT_EQ(seqan3::detail::extract_eighth<3>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 3]);
        EXPECT_EQ(seqan3::detail::extract_eighth<4>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 4]);
        EXPECT_EQ(seqan3::detail::extract_eighth<5>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 5]);
        EXPECT_EQ(seqan3::detail::extract_eighth<6>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 6]);
        EXPECT_EQ(seqan3::detail::extract_eighth<7>(vec)[idx], vec[idx + TestFixture::simd_length / 8 * 7]);
    }
}

//-----------------------------------------------------------------------------
// Algorithm upcast
//-----------------------------------------------------------------------------

template <typename t>
class simd_algorithm_upcast : public ::testing::Test
{
public:

    static constexpr auto target_list_signed()
    {
        if constexpr (sizeof(t) == 1)
            return type_list_signed_8{};
        else if constexpr (sizeof(t) == 2)
            return type_list_signed_16{};
        else
            return type_list_signed_32{};
    }

    static constexpr auto target_list_unsigned()
    {
        if constexpr (sizeof(t) == 1)
            return type_list_unsigned_8{};
        else if constexpr (sizeof(t) == 2)
            return type_list_unsigned_16{};
        else
            return type_list_unsigned_32{};
    }

    using target_list_signed_t = decltype(target_list_signed());
    using target_list_unsigned_t = decltype(target_list_unsigned());

    using type_list_signed_8  = seqan3::type_list<int8_t, int16_t, int32_t, int64_t>;
    using type_list_signed_16 = seqan3::type_list<int16_t, int32_t, int64_t>;
    using type_list_signed_32 = seqan3::type_list<int32_t, int64_t>;

    using type_list_unsigned_8  = seqan3::type_list<uint8_t, uint16_t, uint32_t, uint64_t>;
    using type_list_unsigned_16 = seqan3::type_list<uint16_t, uint32_t, uint64_t>;
    using type_list_unsigned_32 = seqan3::type_list<uint32_t, uint64_t>;
};

using upcast_test_types = ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t>;
TYPED_TEST_SUITE(simd_algorithm_upcast, upcast_test_types, );

TYPED_TEST(simd_algorithm_upcast, signed)
{
    using list = typename TestFixture::target_list_signed_t;

    seqan3::detail::for_each<list>([this] (auto type_id)
    {
        using target_type = typename decltype(type_id)::type;
        using src_simd_t = seqan3::simd::simd_type_t<TypeParam>;
        using target_simd_t = seqan3::simd::simd_type_t<target_type>;

        src_simd_t s = seqan3::simd::fill<src_simd_t>(-10);
        target_simd_t t = seqan3::simd::upcast<target_simd_t>(s);

        // Need to compare elementwise since the simd types are not equal and not comparable with `SIMD_EQ`.
        for (size_t i  = 0; i < seqan3::simd::simd_traits<target_simd_t>::length; ++i)
            EXPECT_EQ(t[i], static_cast<target_type>(static_cast<TypeParam>(-10)));
    });
}

TYPED_TEST(simd_algorithm_upcast, unsigned)
{
    using list = typename TestFixture::target_list_unsigned_t;

    seqan3::detail::for_each<list>([this] (auto type_id)
    {
        using target_type = typename decltype(type_id)::type;
        using src_simd_t = seqan3::simd::simd_type_t<TypeParam>;
        using target_simd_t = seqan3::simd::simd_type_t<target_type>;

        src_simd_t s = seqan3::simd::fill<src_simd_t>(-10);
        target_simd_t t = seqan3::simd::upcast<target_simd_t>(s);

        // Need to compare elementwise since the simd types are not equal and not comparable with `SIMD_EQ`.
        for (size_t i  = 0; i < seqan3::simd::simd_traits<target_simd_t>::length; ++i)
            EXPECT_EQ(t[i], static_cast<target_type>(static_cast<TypeParam>(-10)));
    });
}
