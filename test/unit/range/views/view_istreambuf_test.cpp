// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/istreambuf.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/ranges>
#include <seqan3/test/tmp_filename.hpp>

#include "../iterator_test_template.hpp"

using namespace seqan3;

using iterator_type = decltype(views::istreambuf(std::declval<std::istringstream &>()).begin());

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    std::string expected_range{"ACGTATATATAT ATATAT TTA \n AUAUAA"};
    std::istringstream is{expected_range};

    decltype(views::istreambuf(is)) test_range = views::istreambuf(is);
};

using test_type = ::testing::Types<iterator_type>;

INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

TEST(view_istreambuf, basic)
{
    std::string data{"ACGTATATATAT ATATAT TTA \n AUAUAA"};
    std::istringstream is{data};

    // construct from istream:
    auto v1 = views::istreambuf(is);
    size_t i = 0;
    for (auto c : v1)
        EXPECT_EQ(c, data[i++]);
    EXPECT_EQ(i, data.size());

    // construct from streambuf
    is.clear();
    is.seekg(0, std::ios::beg);
    auto v2 = views::istreambuf(*is.rdbuf());
    i = 0;
    for (auto c : v2)
        EXPECT_EQ(c, data[i++]);
    EXPECT_EQ(i, data.size());

    // combinability
    is.clear();
    is.seekg(0, std::ios::beg);
    auto v3 = views::istreambuf(is) | views::char_to<dna5> | views::complement;
    std::vector<dna5> comp{"TGCATATATATANTATATANAATNNNTATATT"_dna5};
    EXPECT_TRUE(std::ranges::equal(v3, comp));

    // combinability 2
    is.clear();
    is.seekg(0, std::ios::beg);
    auto v4 = views::istreambuf(is) | views::take_until(is_space);
    std::string out2 = v4 | views::to<std::string>;
    std::string comp2 = "ACGTATATATAT";
    EXPECT_TRUE(std::ranges::equal(out2, comp2));
}

TEST(view_istreambuf, concepts)
{
    std::string data{""};
    std::istringstream is{data};
    auto v1 = views::istreambuf(is);

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), char>));
}

TEST(view_istreambuf, big_file_stram)
{
    test::tmp_filename file_name{"istream_storage"};

    {
        std::ofstream os{file_name.get_path()};
        for (size_t idx = 0; idx < 11000 ; ++idx)
            os << "halloballo\n";
    }

    std::ifstream istream{file_name.get_path()};
    auto v = views::istreambuf(istream);
    while (v.begin() != v.end())
    {
        EXPECT_TRUE(std::ranges::equal(v | views::take_until_or_throw_and_consume(is_char<'\n'>),
                                       std::string_view{"halloballo"}));
    }

}
