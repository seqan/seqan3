// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/concat.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation.hpp>

using seqan3::operator""_aa27;
using seqan3::operator""_dna4;
using seqan3::operator""_dna15;

TEST(translate_triplets, dna4)
{
    seqan3::dna4 n1{'C'_dna4};
    seqan3::dna4 n2{'T'_dna4};
    seqan3::dna4 n3{'A'_dna4};
    seqan3::aa27 c{'L'_aa27};

    seqan3::aa27 t1{seqan3::translate_triplet<seqan3::genetic_code::CANONICAL>(n1, n2, n3)};

    EXPECT_EQ(t1, c);
}

TEST(translate_triplets, dna15)
{
    seqan3::dna15 n1{'C'_dna15};
    seqan3::dna15 n2{'T'_dna15};
    seqan3::dna15 n3{'A'_dna15};
    seqan3::aa27 c{'L'_aa27};

    seqan3::aa27 t1{seqan3::translate_triplet<seqan3::genetic_code::CANONICAL, seqan3::dna15>(n1, n2, n3)};

    EXPECT_EQ(t1, c);
}

TEST(translate_triplets, random_access_range)
{
    seqan3::dna15 n1{'C'_dna15};
    seqan3::dna15 n2{'T'_dna15};
    seqan3::dna15 n3{'A'_dna15};
    seqan3::aa27 c{'L'_aa27};

    auto range_triplet = ranges::view::concat(std::views::single(n1), std::views::single(n2),
                                              std::views::single(n3));
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    seqan3::aa27 t2{seqan3::translate_triplet(range_triplet)};

    EXPECT_TRUE(std::ranges::random_access_range<decltype(range_triplet)>);
    EXPECT_EQ(t2, c);
}

TEST(translate_triplets, tuple)
{
    seqan3::dna15 n1{'C'_dna15};
    seqan3::dna15 n2{'T'_dna15};
    seqan3::dna15 n3{'A'_dna15};
    seqan3::aa27 c{'L'_aa27};

    std::tuple tuple_triplet{n1, n2, n3};
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    seqan3::aa27 t3{seqan3::translate_triplet(tuple_triplet)};

    EXPECT_EQ(t3, c);
}
