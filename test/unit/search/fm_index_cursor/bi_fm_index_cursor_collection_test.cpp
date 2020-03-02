// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "bi_fm_index_cursor_collection_test_template.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/to.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

// dna4 - bi_fm_index
template <typename sdsl_type>
struct bi_fm_index_cursor_collection_test<seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4,
                                                                                         seqan3::text_layout::collection,
                                                                                         sdsl_type>>> :
    public ::testing::Test
{
    const seqan3::dna4_vector text{"ACGGTAGGACGTAGC"_dna4};
    const seqan3::dna4_vector text1{"AACGATCGGA"_dna4};
    const seqan3::dna4_vector text2{"TGCTACGATCC"_dna4};
    const seqan3::dna4_vector text3 = seqan3::views::slice(text, 0, 11)
                              | seqan3::views::to<seqan3::dna4_vector>; // "ACGGTAGGACG"
    const seqan3::dna4_vector text4 = seqan3::views::slice(text, 0, 14)
                              | seqan3::views::to<seqan3::dna4_vector>; // "ACGGTAGGACGTAG"

    const std::vector<seqan3::dna4_vector> text_col1{text1, text1};
    const std::vector<seqan3::dna4_vector> text_col2{text3, text2};
    const std::vector<seqan3::dna4_vector> text_col3{text4, text2};
    const std::vector<seqan3::dna4_vector> text_col4{text, text2};

    const std::vector<seqan3::dna4_vector> rev_text1 = text_col1 | seqan3::views::deep{std::views::reverse}
                                                           | seqan3::views::deep{seqan3::views::persist}
                                                           | seqan3::views::to<const std::vector<seqan3::dna4_vector>>;
    const std::vector<seqan3::dna4_vector> rev_text2 = text_col4 | seqan3::views::deep{std::views::reverse}
                                                           | std::views::reverse
                                                           | seqan3::views::to<const std::vector<seqan3::dna4_vector>>;

    const seqan3::dna4_vector pattern1{"CAG"_dna4};
    const seqan3::dna4_vector pattern2{"TT"_dna4};
    const seqan3::dna4_vector pattern3{"GATGC"_dna4};
    const seqan3::dna4_vector pattern4{"GATGG"_dna4};
};

// dna5 - bi_fm_index
template <typename sdsl_type>
struct bi_fm_index_cursor_collection_test<seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna5,
                                                                                         seqan3::text_layout::collection,
                                                                                         sdsl_type>>> :
    public ::testing::Test
{
    const seqan3::dna5_vector text{"ACGGTAGGACGTAGC"_dna5};
    const seqan3::dna5_vector text1{"AACGATCGGA"_dna5};
    const seqan3::dna5_vector text2{"TGCTACGATCC"_dna5};
    const seqan3::dna5_vector text3 = seqan3::views::slice(text, 0, 11)
                              | seqan3::views::to<seqan3::dna5_vector>; // "ACGGTAGGACG"
    const seqan3::dna5_vector text4 = seqan3::views::slice(text, 0, 14)
                              | seqan3::views::to<seqan3::dna5_vector>; // "ACGGTAGGACGTAG"

    const std::vector<seqan3::dna5_vector> text_col1{text1, text1};
    const std::vector<seqan3::dna5_vector> text_col2{text3, text2};
    const std::vector<seqan3::dna5_vector> text_col3{text4, text2};
    const std::vector<seqan3::dna5_vector> text_col4{text, text2};

    const std::vector<seqan3::dna5_vector> rev_text1 = text_col1 | seqan3::views::deep{std::views::reverse}
                                                           | seqan3::views::deep{seqan3::views::persist}
                                                           | seqan3::views::to<const std::vector<seqan3::dna5_vector>>;
    const std::vector<seqan3::dna5_vector> rev_text2 = text_col4 | seqan3::views::deep{std::views::reverse}
                                                           | std::views::reverse
                                                           | seqan3::views::to<const std::vector<seqan3::dna5_vector>>;

    const seqan3::dna5_vector pattern1{"CAG"_dna5};
    const seqan3::dna5_vector pattern2{"TT"_dna5};
    const seqan3::dna5_vector pattern3{"GATGC"_dna5};
    const seqan3::dna5_vector pattern4{"GATGG"_dna5};
};

using it_t1 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, bi_fm_index_cursor_collection_test, it_t1, );

using it_t2 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, bi_fm_index_cursor_collection_test, it_t2, );
