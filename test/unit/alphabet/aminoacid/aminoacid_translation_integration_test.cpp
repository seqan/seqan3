// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <range/v3/view/concat.hpp>
#include <range/v3/view/single.hpp>

#include <seqan3/alphabet/aminoacid/all.hpp>

using namespace seqan3;

TEST(translation, translate_triplets)
{
    dna15 n1{'C'_dna15};
    dna15 n2{'T'_dna15};
    dna15 n3{'A'_dna15};
    aa27 c{aa27::L};

    // Nucleotide interface
    aa27 t1{translate_triplet<genetic_code::CANONICAL, dna15>(n1, n2, n3)};

    EXPECT_EQ(t1, c);

    // Range interface
    auto range_triplet = ranges::view::concat(ranges::view::single(n1), ranges::view::single(n2),
                                              ranges::view::single(n3));
    aa27 t2{translate_triplet(range_triplet)};

    EXPECT_EQ(t2, c);

    // Tuple interface
    std::tuple tuple_triplet{n1, n2, n3};
    aa27 t3{translate_triplet(tuple_triplet)};

    EXPECT_EQ(t3, c);
}
