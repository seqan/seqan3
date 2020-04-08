// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_collection_test_template.hpp"
#include "fm_index_test_template.hpp"

using t1 = std::pair<seqan3::bi_fm_index<unsigned char, seqan3::text_layout::single>,
                     std::vector<unsigned char>>;
INSTANTIATE_TYPED_TEST_SUITE_P(char, fm_index_test, t1, );
using t2 = std::pair<seqan3::bi_fm_index<unsigned char, seqan3::text_layout::collection>,
                     std::vector<std::vector<unsigned char>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(char_collection, fm_index_collection_test, t2, );

TEST(char, throw_on_reserved_char)
{
    using bi_fm_index_t = seqan3::bi_fm_index<unsigned char, seqan3::text_layout::single>;

    unsigned char c = 255;
    std::vector<unsigned char> text{'a', 'u', ',', c, '0'};

    EXPECT_THROW(bi_fm_index_t bi_index{text}, std::out_of_range);
}

TEST(char_collection, throw_on_reserved_char)
{
    using bi_fm_index_t = seqan3::bi_fm_index<unsigned char, seqan3::text_layout::collection>;

    {
        unsigned char c = 255;
        std::vector<std::vector<unsigned char>> text{{'a', 'b'}, {'a', 'u', ',', c, '0'}};

        EXPECT_THROW(bi_fm_index_t bi_index{text}, std::out_of_range);
    }
    {
        unsigned char c = 254;
        std::vector<std::vector<unsigned char>> text{{'a', 'b'}, {'a', 'u', ',', c, '0'}};

        EXPECT_THROW(bi_fm_index_t bi_index{text}, std::out_of_range);
    }
}
