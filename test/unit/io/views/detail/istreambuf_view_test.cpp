// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <ranges>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include "../../../range/iterator_test_template.hpp"

using seqan3::operator""_dna5;

using iterator_type = decltype(seqan3::detail::istreambuf(std::declval<std::istringstream &>()).begin());

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    std::string expected_range{"ACGTATATATAT ATATAT TTA \n AUAUAA"};
    std::istringstream is{expected_range};

    decltype(seqan3::detail::istreambuf(is)) test_range = seqan3::detail::istreambuf(is);
};

using test_type = ::testing::Types<iterator_type>;

INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

TEST(view_istreambuf, basic)
{
    using namespace std::literals;

    std::string data{"ACGTATATATAT ATATAT TTA \n AUAUAA"};
    std::istringstream is{data};

    // construct from istream:
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(is), data)

    // construct from streambuf
    is.clear();
    is.seekg(0, std::ios::beg);
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(*is.rdbuf()), data)

    // combinability
    is.clear();
    is.seekg(0, std::ios::beg);
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(is) | seqan3::views::char_to<seqan3::dna5> | seqan3::views::complement,
                    "TGCATATATATANTATATANAATNNNTATATT"_dna5);

    // combinability 2
    is.clear();
    is.seekg(0, std::ios::beg);
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(is) | seqan3::detail::take_until(seqan3::is_space), "ACGTATATATAT"sv);
}

TEST(view_istreambuf, concepts)
{
    std::string data{""};
    std::istringstream is{data};
    auto v1 = seqan3::detail::istreambuf(is);

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), char>));
}

TEST(view_istreambuf, big_file_stram)
{
    using namespace std::literals;

    seqan3::test::tmp_directory tmp_dir{};
    auto file_name = tmp_dir.path() / "istream_storage";

    {
        std::ofstream os{file_name};
        for (size_t idx = 0; idx < 11000; ++idx)
            os << "halloballo\n";
    }

    std::ifstream istream{file_name};
    auto v = seqan3::detail::istreambuf(istream);
    while (v.begin() != v.end())
    {
        EXPECT_RANGE_EQ(v | seqan3::detail::take_until_or_throw_and_consume(seqan3::is_char<'\n'>), "halloballo"sv);
    }
}
