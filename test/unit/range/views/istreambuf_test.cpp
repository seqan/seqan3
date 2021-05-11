// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/io/detail/istreambuf_view.hpp>
#include <seqan3/io/detail/take_until_view.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include "../iterator_test_template.hpp"

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
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(is),
                    data)

    // construct from streambuf
    is.clear();
    is.seekg(0, std::ios::beg);
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(*is.rdbuf()),
                    data)

    // combinability
    is.clear();
    is.seekg(0, std::ios::beg);
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(is) | seqan3::views::char_to<seqan3::dna5> | seqan3::views::complement,
                    "TGCATATATATANTATATANAATNNNTATATT"_dna5);

    // combinability 2
    is.clear();
    is.seekg(0, std::ios::beg);
    EXPECT_RANGE_EQ(seqan3::detail::istreambuf(is) | seqan3::detail::take_until(seqan3::is_space),
                    "ACGTATATATAT"sv);
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

    seqan3::test::tmp_filename file_name{"istream_storage"};

    {
        std::ofstream os{file_name.get_path()};
        for (size_t idx = 0; idx < 11000 ; ++idx)
            os << "halloballo\n";
    }

    std::ifstream istream{file_name.get_path()};
    auto v = seqan3::detail::istreambuf(istream);
    while (v.begin() != v.end())
    {
        EXPECT_RANGE_EQ(v | seqan3::detail::take_until_or_throw_and_consume(seqan3::is_char<'\n'>),
                        "halloballo"sv);
    }

}
