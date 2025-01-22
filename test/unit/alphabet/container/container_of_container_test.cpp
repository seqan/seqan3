// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/cereal.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>

using seqan3::operator""_dna4;

template <typename T>
class container_of_container : public ::testing::Test
{};

using container_of_container_types =
    ::testing::Types<std::vector<std::vector<seqan3::dna4>>,
                     seqan3::concatenated_sequences<std::vector<seqan3::dna4>>,
                     seqan3::concatenated_sequences<seqan3::bitpacked_sequence<seqan3::dna4>>>;

TYPED_TEST_SUITE(container_of_container, container_of_container_types, );

TYPED_TEST(container_of_container, concepts)
{
    EXPECT_TRUE(seqan3::container<TypeParam>);
}

TYPED_TEST(container_of_container, construction)
{
    TypeParam t1;
    TypeParam t2{};
    EXPECT_EQ(t1, t2);

    // initializer list
    TypeParam t3{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam t4 = {"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    EXPECT_EQ(t3, t4);

    // n * value
    TypeParam t5{2, "ACGT"_dna4};
    // from another TypeParam's sub-range
    TypeParam t6{t3.begin(), t3.begin() + 2};
    EXPECT_EQ(t5, t6);

    std::vector<std::vector<seqan3::dna4>> other_vector{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    // direct from another container-of-container
    TypeParam t7{other_vector};
    // from another container-of-container's sub-range
    TypeParam t8{other_vector.cbegin(), other_vector.cend()};
    EXPECT_EQ(t3, t7);
    EXPECT_EQ(t7, t8);
}

TYPED_TEST(container_of_container, assign)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam t2{"ACGT"_dna4, "ACGT"_dna4};
    std::vector<std::vector<seqan3::dna4>> other_vector{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // n * value
    TypeParam t3;
    t3.assign(2, "ACGT"_dna4);
    EXPECT_EQ(t3, t2);

    // from another container-of-container's sub-range
    TypeParam t4;
    t4.assign(other_vector.cbegin(), other_vector.cend());
    EXPECT_EQ(t4, t1);

    // initializer list
    TypeParam t5, t6;
    t5.assign({"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4});
    t6 = {"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    EXPECT_EQ(t5, t1);
    EXPECT_EQ(t6, t1);

    // direct from another container-of-container
    if constexpr (std::is_same_v<TypeParam, seqan3::concatenated_sequences<std::vector<seqan3::dna4>>>)
    {
        TypeParam t7, t8;
        t7.assign(other_vector);
        t8 = other_vector;
        EXPECT_EQ(t7, t1);
        EXPECT_EQ(t8, t1);
    }
}

TYPED_TEST(container_of_container, iterators)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam const t2{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // begin
    EXPECT_RANGE_EQ(*t1.begin(), "ACGT"_dna4);
    EXPECT_RANGE_EQ(*t1.cbegin(), "ACGT"_dna4);
    EXPECT_RANGE_EQ(*t2.begin(), "ACGT"_dna4);
    EXPECT_RANGE_EQ(*t2.cbegin(), "ACGT"_dna4);

    // end and arithmetic
    EXPECT_RANGE_EQ(*(t1.end() - 1), "GAGGA"_dna4);
    EXPECT_RANGE_EQ(*(t1.cend() - 1), "GAGGA"_dna4);
    EXPECT_RANGE_EQ(*(t2.end() - 1), "GAGGA"_dna4);
    EXPECT_RANGE_EQ(*(t2.cend() - 1), "GAGGA"_dna4);

    // convertibility between const and non-const
    EXPECT_TRUE(t1.cend() == t1.end());

    // writability
    (*t1.begin())[0] = 'T'_dna4;
    EXPECT_RANGE_EQ(*t1.begin(), "TCGT"_dna4);
}

TYPED_TEST(container_of_container, element_access)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam const t2{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // at
    EXPECT_RANGE_EQ(t1.at(0), "ACGT"_dna4);
    EXPECT_RANGE_EQ(t2.at(0), "ACGT"_dna4);
    //TODO once we have throwing assert, check at's ability to throw

    // []
    EXPECT_RANGE_EQ(t1[0], "ACGT"_dna4);
    EXPECT_RANGE_EQ(t2[0], "ACGT"_dna4);

    // front
    EXPECT_RANGE_EQ(t1.front(), "ACGT"_dna4);
    EXPECT_RANGE_EQ(t2.front(), "ACGT"_dna4);

    // back
    EXPECT_RANGE_EQ(t1.back(), "GAGGA"_dna4);
    EXPECT_RANGE_EQ(t2.back(), "GAGGA"_dna4);

    if constexpr (std::is_same_v<TypeParam, seqan3::concatenated_sequences<std::vector<seqan3::dna4>>>)
    {
        using size_type = typename TypeParam::size_type;
        // concat
        EXPECT_RANGE_EQ(t1.concat(), "ACGTACGTGAGGA"_dna4);
        EXPECT_RANGE_EQ(t2.concat(), "ACGTACGTGAGGA"_dna4);

        // data
        EXPECT_EQ(std::get<0>(t1.raw_data()), "ACGTACGTGAGGA"_dna4);
        EXPECT_EQ(std::get<0>(t2.raw_data()), "ACGTACGTGAGGA"_dna4);
        EXPECT_EQ(std::get<1>(t1.raw_data()), (std::vector<size_type>{0, 4, 8, 13}));
        EXPECT_EQ(std::get<1>(t2.raw_data()), (std::vector<size_type>{0, 4, 8, 13}));
    }
}

TYPED_TEST(container_of_container, capacity)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    TypeParam const t2{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // empty
    EXPECT_TRUE(t0.empty());
    EXPECT_FALSE(t1.empty());
    EXPECT_FALSE(t2.empty());

    // size
    EXPECT_EQ(t0.size(), 0u);
    EXPECT_EQ(t1.size(), 3u);
    EXPECT_EQ(t2.size(), 3u);

    // max_size
    EXPECT_GT(t0.max_size(), 1'000'000'000'000u);
    EXPECT_GT(t1.max_size(), 1'000'000'000'000u);
    EXPECT_GT(t2.max_size(), 1'000'000'000'000u);

    // capacity
    EXPECT_GE(t0.capacity(), t0.size());
    EXPECT_GE(t1.capacity(), t1.size());
    EXPECT_GE(t2.capacity(), t2.size());

    // reserve
    EXPECT_LT(t0.capacity(), 1000u);
    t0.reserve(1000);
    EXPECT_GE(t0.capacity(), 1000u);

    // shrink_to_fit
    t1.reserve(1000);
    EXPECT_GT(t1.capacity(), t1.size() * 2);
    t1.shrink_to_fit();
    EXPECT_LE(t1.capacity(), t1.size() * 2);

    if constexpr (std::is_same_v<TypeParam, seqan3::concatenated_sequences<std::vector<seqan3::dna4>>>)
    {
        // size
        EXPECT_EQ(t0.concat_size(), 0u);
        EXPECT_EQ(t1.concat_size(), 13u);
        EXPECT_EQ(t2.concat_size(), 13u);

        // capacity
        EXPECT_GE(t0.concat_capacity(), t0.concat_size());
        EXPECT_GE(t1.concat_capacity(), t1.concat_size());
        EXPECT_GE(t2.concat_capacity(), t2.concat_size());

        // reserve
        EXPECT_LT(t0.concat_capacity(), 1000u);
        t0.concat_reserve(1000);
        EXPECT_GE(t0.concat_capacity(), 1000u);
    }
}

TYPED_TEST(container_of_container, clear)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    t1.clear();
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container_of_container, insert)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // position, value
    t0.insert(t0.cend(), "ACGT"_dna4);
    t0.insert(t0.cend(), "GAGGA"_dna4);
    t0.insert(t0.cbegin() + 1, "ACGT"_dna4);
    EXPECT_EQ(t0, t1);

    // position, n times values
    t0.clear();
    t1 = {"GAGGA"_dna4, "ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t0.insert(t0.cend(), 2, "ACGT"_dna4);
    t0.insert(t0.cend(), 1, "GAGGA"_dna4);
    t0.insert(t0.cbegin(), 1, "GAGGA"_dna4);
    EXPECT_EQ(t0, t1);

    // iterator pair
    t0.clear();
    t1 = {"GAGGA"_dna4, "ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t0.insert(t0.cend(), t1.begin() + 1, t1.begin() + 3);

    t0.insert(t0.cend(), t1.cend() - 1, t1.cend());
    t0.insert(t0.cbegin(), t1.cend() - 1, t1.cend());
    EXPECT_EQ(t0, t1);

    // initializer list
    t0.clear();
    t1 = {"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t0.insert(t0.cend(), {"ACGT"_dna4, "GAGGA"_dna4});
    t0.insert(t0.cbegin() + 1, "ACGT"_dna4);
    EXPECT_EQ(t0, t1);
}

TYPED_TEST(container_of_container, erase)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    // one element
    t1.erase(t1.begin());
    EXPECT_EQ(t1, (TypeParam{"ACGT"_dna4, "GAGGA"_dna4}));

    // range
    t1 = {"GAGGA"_dna4, "ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    t1.erase(t1.begin() + 1, t1.begin() + 3);
    EXPECT_EQ(t1, (TypeParam{"GAGGA"_dna4, "GAGGA"_dna4}));
}

TYPED_TEST(container_of_container, push_pop)
{
    TypeParam t0{};

    // push_back
    t0.push_back("ACGT"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4}));
    t0.push_back("GAGGA"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGA"_dna4}));

    // pop_back
    t0.pop_back();
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4}));
    t0.pop_back();
    EXPECT_EQ(t0, TypeParam{});
}

TYPED_TEST(container_of_container, resize)
{
    TypeParam t0{};

    // enlarge without values
    t0.resize(3);
    EXPECT_EQ(t0, (TypeParam{{}, {}, {}}));

    // enlarge with value
    t0.resize(5, "ACGT"_dna4);
    EXPECT_EQ(t0, (TypeParam{{}, {}, {}, "ACGT"_dna4, "ACGT"_dna4}));

    // shrink with value
    t0.resize(4, "ACGT"_dna4);
    EXPECT_EQ(t0, (TypeParam{{}, {}, {}, "ACGT"_dna4}));

    // shrink without value
    t0.resize(2);
    EXPECT_EQ(t0, (TypeParam{{}, {}}));
}

TYPED_TEST(container_of_container, swap)
{
    TypeParam t0{};
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};

    t0.swap(t1);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4}));
    EXPECT_EQ(t1, TypeParam{});
}

TYPED_TEST(container_of_container, serialisation)
{
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    seqan3::test::do_serialisation(t1);
}

template <typename T>
class concatenated_sequences : public ::testing::Test
{};

using concatenated_sequences_types =
    ::testing::Types<seqan3::concatenated_sequences<std::vector<seqan3::dna4>>,
                     seqan3::concatenated_sequences<seqan3::bitpacked_sequence<seqan3::dna4>>>;

TYPED_TEST_SUITE(concatenated_sequences, concatenated_sequences_types, );

TYPED_TEST(concatenated_sequences, last_push_back)
{
    TypeParam t0{};

    t0.push_back("ACGT"_dna4);
    t0.push_back("GAGGA"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGA"_dna4}));

    t0.last_push_back('C'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGAC"_dna4}));
    t0.last_push_back('G'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGACG"_dna4}));
    t0.last_push_back('T'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGACGT"_dna4}));

    t0.push_back("ACGT"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGACGT"_dna4, "ACGT"_dna4}));
    t0.last_push_back('C'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGACGT"_dna4, "ACGTC"_dna4}));
    t0.last_push_back('G'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGACGT"_dna4, "ACGTCG"_dna4}));
    t0.last_push_back('T'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "GAGGACGT"_dna4, "ACGTCGT"_dna4}));
}

TYPED_TEST(concatenated_sequences, push_back_empty)
{
    TypeParam t0{};

    t0.push_back("ACGT"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4}));

    t0.push_back();
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, ""_dna4}));
    t0.last_push_back('C'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "C"_dna4}));
    t0.last_push_back('G'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CG"_dna4}));
    t0.last_push_back('T'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4}));

    t0.push_back();
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4, ""_dna4}));
    t0.last_push_back('C'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4, "C"_dna4}));
    t0.last_push_back('G'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4, "CG"_dna4}));
    t0.last_push_back('T'_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4, "CGT"_dna4}));
}

TYPED_TEST(concatenated_sequences, last_append)
{
    TypeParam t0{};

    t0.push_back("ACGT"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4}));

    t0.push_back();
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, ""_dna4}));
    t0.last_append("C"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "C"_dna4}));
    t0.last_append("GT"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4}));

    t0.push_back();
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4, ""_dna4}));
    t0.last_append("C"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4, "C"_dna4}));
    t0.last_append("GT"_dna4);
    EXPECT_EQ(t0, (TypeParam{"ACGT"_dna4, "CGT"_dna4, "CGT"_dna4}));
}

TEST(concatenated_sequences_, associated_types)
{
    {
        using t = seqan3::concatenated_sequences<std::vector<int>>;

        EXPECT_TRUE((std::same_as<std::ranges::range_value_t<t>, std::span<int>>));
        EXPECT_TRUE((std::same_as<std::ranges::range_reference_t<t>, std::span<int>>));

        EXPECT_TRUE((std::same_as<std::ranges::range_value_t<t const>, std::span<int>>));
        EXPECT_TRUE((std::same_as<std::ranges::range_reference_t<t const>, std::span<int const>>));
    }

    {
        using t = seqan3::concatenated_sequences<std::string>;

        EXPECT_TRUE((std::same_as<std::ranges::range_value_t<t>, std::span<char>>));
        EXPECT_TRUE((std::same_as<std::ranges::range_reference_t<t>, std::span<char>>));

        EXPECT_TRUE((std::same_as<std::ranges::range_value_t<t const>, std::span<char>>));
        EXPECT_TRUE((std::same_as<std::ranges::range_reference_t<t const>, std::string_view>));
    }

    {
        using t = seqan3::concatenated_sequences<seqan3::bitpacked_sequence<seqan3::dna4>>;
        using it_t = std::ranges::iterator_t<seqan3::bitpacked_sequence<seqan3::dna4>>;
        using sen_t = std::ranges::sentinel_t<seqan3::bitpacked_sequence<seqan3::dna4>>;

        using c_it_t = std::ranges::iterator_t<seqan3::bitpacked_sequence<seqan3::dna4> const>;
        using c_sen_t = std::ranges::sentinel_t<seqan3::bitpacked_sequence<seqan3::dna4> const>;

        EXPECT_TRUE((std::same_as<std::ranges::range_value_t<t>, std::ranges::subrange<it_t, sen_t>>));
        EXPECT_TRUE((std::same_as<std::ranges::range_reference_t<t>, std::ranges::subrange<it_t, sen_t>>));

        EXPECT_TRUE((std::same_as<std::ranges::range_value_t<t const>, std::ranges::subrange<it_t, sen_t>>));
        EXPECT_TRUE((std::same_as<std::ranges::range_reference_t<t const>, std::ranges::subrange<c_it_t, c_sen_t>>));
    }
}
