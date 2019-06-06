// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/range/view/istreambuf.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/ranges>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

TEST(view_istreambuf, basic)
{
    std::string data{"ACGTATATATAT ATATAT TTA \n AUAUAA"};
    std::istringstream is{data};

    // construct from istream:
    auto v1 = view::istreambuf(is);
    size_t i = 0;
    for (auto c : v1)
        EXPECT_EQ(c, data[i++]);
    EXPECT_EQ(i, data.size());

    // construct from streambuf
    is.clear();
    is.seekg(0, std::ios::beg);
    auto v2 = view::istreambuf(*is.rdbuf());
    i = 0;
    for (auto c : v2)
        EXPECT_EQ(c, data[i++]);
    EXPECT_EQ(i, data.size());

    // combinability
    is.clear();
    is.seekg(0, std::ios::beg);
    auto v3 = view::istreambuf(is) | view::char_to<dna5> | view::complement;
    std::vector<dna5> comp{"TGCATATATATANTATATANAATNNNTATATT"_dna5};
    EXPECT_TRUE(std::ranges::equal(v3, comp));

    // combinability 2
    is.clear();
    is.seekg(0, std::ios::beg);
    auto v4 = view::istreambuf(is) | view::take_until(is_space);
    std::string out2 = v4;
    std::string comp2 = "ACGTATATATAT";
    EXPECT_TRUE(std::ranges::equal(out2, comp2));
}

TEST(view_istreambuf, concepts)
{
    std::string data{""};
    std::istringstream is{data};
    auto v1 = view::istreambuf(is);

    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_FALSE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_FALSE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_FALSE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_FALSE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(ConstIterableRange<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), char>));
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
    auto v = view::istreambuf(istream);
    while (v.begin() != v.end())
    {
        EXPECT_TRUE(std::ranges::equal(v | view::take_until_or_throw_and_consume(is_char<'\n'>),
                                       std::string_view{"halloballo"}));
    }

}
