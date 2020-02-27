// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_cursor_collection_test_template.hpp"

#include <climits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

// dna4 - fm_index
template <typename sdsl_type>
struct fm_index_cursor_collection_test<seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4,
                                                                                seqan3::text_layout::collection,
                                                                                sdsl_type>>> : public ::testing::Test
{
    using alphabet_type = seqan3::dna4;

    const seqan3::dna4_vector text1{"ACGACG"_dna4};
    const seqan3::dna4_vector text2{"ACGAACGC"_dna4};
    const seqan3::dna4_vector text3{"CGTCGT"_dna4};
    const seqan3::dna4_vector text4{"ATATAT"_dna4};
    const seqan3::dna4_vector text5{"TGCGATCGA"_dna4};
    const seqan3::dna4_vector text6{"TACGATCGA"_dna4};
    const seqan3::dna4_vector text7{"ACGTACGT"_dna4};
    const seqan3::dna4_vector text8{"TGCGATACGA"_dna4};

    seqan3::dna4_vector empty_text{""_dna4};
    const size_t max_rank = seqan3::dna4::alphabet_size -1; // T

    const std::vector<seqan3::dna4_vector> text_col1{text1, text1};
    const std::vector<seqan3::dna4_vector> text_col2{text1, text5};
    const std::vector<seqan3::dna4_vector> text_col3{text1, empty_text, empty_text, text5};
    const std::vector<seqan3::dna4_vector> text_col4{text2, text6};
    const std::vector<seqan3::dna4_vector> text_col5{text2, text5};
    const std::vector<seqan3::dna4_vector> text_col6{text3, text3};
    const std::vector<seqan3::dna4_vector> text_col7{text4, text4};
    const std::vector<seqan3::dna4_vector> text_col8{text7, text8};
};

// dna4 - bi_fm_index
template <typename alphabet_type, typename sdsl_type>
struct fm_index_cursor_collection_test<seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<alphabet_type,
                                                                                      seqan3::text_layout::collection,
                                                                                      sdsl_type>>> :
    public fm_index_cursor_collection_test<seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4,
                                                                                    seqan3::text_layout::collection,
                                                                                    sdsl_type>>>
{};

// dna5 - fm_index
template <typename sdsl_type>
struct fm_index_cursor_collection_test<seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna5,
                                                                                seqan3::text_layout::collection,
                                                                                sdsl_type>>> : public ::testing::Test
{
    using alphabet_type = seqan3::dna5;

    const seqan3::dna5_vector text1{"ACGACG"_dna5};
    const seqan3::dna5_vector text2{"ACGAACGC"_dna5};
    const seqan3::dna5_vector text3{"CGTCGT"_dna5};
    const seqan3::dna5_vector text4{"ATATAT"_dna5};
    const seqan3::dna5_vector text5{"TGCGATCGA"_dna5};
    const seqan3::dna5_vector text6{"TACGATCGA"_dna5};
    const seqan3::dna5_vector text7{"ACGTACGT"_dna5};
    const seqan3::dna5_vector text8{"TGCGATACGA"_dna5};

    seqan3::dna5_vector empty_text{""_dna5};
    const size_t max_rank = seqan3::dna5::alphabet_size -1; // T

    const std::vector<seqan3::dna5_vector> text_col1{text1, text1};
    const std::vector<seqan3::dna5_vector> text_col2{text1, text5};
    const std::vector<seqan3::dna5_vector> text_col3{text1, empty_text, empty_text, text5};
    const std::vector<seqan3::dna5_vector> text_col4{text2, text6};
    const std::vector<seqan3::dna5_vector> text_col5{text2, text5};
    const std::vector<seqan3::dna5_vector> text_col6{text3, text3};
    const std::vector<seqan3::dna5_vector> text_col7{text4, text4};
    const std::vector<seqan3::dna5_vector> text_col8{text7, text8};
};

// char - fm_index
template <typename sdsl_type>
struct fm_index_cursor_collection_test<seqan3::fm_index_cursor<seqan3::fm_index<char,
                                                                                seqan3::text_layout::collection,
                                                                                sdsl_type>>> : public ::testing::Test
{
    using alphabet_type = char;

    const std::string text1{"ACGACG"};
    const std::string text2{"ACGAACGC"};
    const std::string text3{"CGTCGT"};
    const std::string text4{"ATATAT"};
    const std::string text5{"TGCGATCGA"};
    const std::string text6{"TACGATCGA"};
    const std::string text7{"ACGTACGT"};
    const std::string text8{"TGCGATACGA"};

    const std::string empty_text{""};
    const size_t max_rank = CHAR_MAX;

     const std::vector<std::string> text_col1{text1, text1};
     const std::vector<std::string> text_col2{text1, text5};
     const std::vector<std::string> text_col3{text1, empty_text, empty_text, text5};
     const std::vector<std::string> text_col4{text2, text6};
     const std::vector<std::string> text_col5{text2, text5};
     const std::vector<std::string> text_col6{text3, text3};
     const std::vector<std::string> text_col7{text4, text4};
     const std::vector<std::string> text_col8{text6, text7};
};

// dna4
using it_t1 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(default_traits, fm_index_cursor_collection_test, it_t1, );

using it_t2 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4,
                                                       seqan3::text_layout::collection,
                                                       sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(byte_alphabet_traits, fm_index_cursor_collection_test, it_t2, );

using it_t3 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_default_traits, fm_index_cursor_collection_test, it_t3, );

using it_t4 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4,
                                                             seqan3::text_layout::collection,
                                                             sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_byte_alphabet_traits, fm_index_cursor_collection_test, it_t4, );

// dna5
using it_t5 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna5_default_traits, fm_index_cursor_collection_test, it_t5, );

// char
using it_t6 = seqan3::fm_index_cursor<seqan3::fm_index<char, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(char_default_traits, fm_index_cursor_collection_test, it_t6, );
