// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "bi_fm_index_cursor_test_template.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/to.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

// dna4 - bi_fm_index
template <typename sdsl_type>
struct bi_fm_index_cursor_test<seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4,
                                                                              seqan3::text_layout::single,
                                                                              sdsl_type>>> : public ::testing::Test
{
    const seqan3::dna4_vector text{"ACGGTAGGACGTAGC"_dna4};
    const seqan3::dna4_vector text1{"AACGATCGGA"_dna4};

    std::vector<seqan3::dna4> rev_text1 = text | std::views::reverse | seqan3::views::to<seqan3::dna4_vector>;
    std::vector<seqan3::dna4> rev_text2 = std::views::reverse(text1) | seqan3::views::to<seqan3::dna4_vector>;

    const seqan3::dna4_vector pattern1{"CAG"_dna4};
    const seqan3::dna4_vector pattern2{"TT"_dna4};
    const seqan3::dna4_vector pattern3{"GATGC"_dna4};
    const seqan3::dna4_vector pattern4{"GATGG"_dna4};
};

// dna5 - bi_fm_index
template <typename sdsl_type>
struct bi_fm_index_cursor_test<seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna5,
                                                                              seqan3::text_layout::single,
                                                                              sdsl_type>>> : public ::testing::Test
{
    const seqan3::dna5_vector text{"ACGGTAGGACGTAGC"_dna5};
    const seqan3::dna5_vector text1{"AACGATCGGA"_dna5};

    std::vector<seqan3::dna5> rev_text1 = text | std::views::reverse | seqan3::views::to<seqan3::dna5_vector>;
    std::vector<seqan3::dna5> rev_text2 = std::views::reverse(text1) | seqan3::views::to<seqan3::dna5_vector>;

    const seqan3::dna5_vector pattern1{"CAG"_dna5};
    const seqan3::dna5_vector pattern2{"TT"_dna5};
    const seqan3::dna5_vector pattern3{"GATGC"_dna5};
    const seqan3::dna5_vector pattern4{"GATGG"_dna5};
};

// char - bi_fm_index
template <typename sdsl_type>
struct bi_fm_index_cursor_test<seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<char,
                                                                              seqan3::text_layout::single,
                                                                              sdsl_type>>> : public ::testing::Test
{
    const std::string text{"ACGGTAGGACGTAGC"};
    const std::string text1{"AACGATCGGA"};

    std::string rev_text1 = text | std::views::reverse | seqan3::views::to<std::string>;
    std::string rev_text2 = std::views::reverse(text1) | seqan3::views::to<std::string>;

    const std::string pattern1{"CAG"};
    const std::string pattern2{"TT"};
    const std::string pattern3{"GATGC"};
    const std::string pattern4{"GATGG"};
};

// dna4
using it_t1 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, bi_fm_index_cursor_test, it_t1, );

// dna5
using it_t2 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, bi_fm_index_cursor_test, it_t2, );

// char
using it_t3 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<char, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(char, bi_fm_index_cursor_test, it_t3, );
