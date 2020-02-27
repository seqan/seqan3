// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_cursor_test_template.hpp"

#include <climits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

// dna4 - fm_index
template <typename sdsl_type>
struct fm_index_cursor_test<seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4,
                                                                    seqan3::text_layout::single,
                                                                    sdsl_type>>> : public ::testing::Test
{
    using alphabet_type = seqan3::dna4;

    const seqan3::dna4_vector text1{"ACGACG"_dna4};
    const seqan3::dna4_vector text2{"ACGAACGC"_dna4};
    const seqan3::dna4_vector text3{"CGTCGT"_dna4};
    const seqan3::dna4_vector text4{"ATATAT"_dna4};
    const seqan3::dna4_vector empty_text{""_dna4};
    const size_t max_rank = seqan3::dna4::alphabet_size -1; // T
};

// dna4 - bi_fm_index
template <typename sdsl_type>
struct fm_index_cursor_test<seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4,
                                                                           seqan3::text_layout::single,
                                                                           sdsl_type>>> :
    public fm_index_cursor_test<seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4,
                                                                         seqan3::text_layout::single,
                                                                         sdsl_type>>>
{};

// dna5 - fm_index
template <typename sdsl_type>
struct fm_index_cursor_test<seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna5,
                                                                    seqan3::text_layout::single,
                                                                    sdsl_type>>> : public ::testing::Test
{
    using alphabet_type = seqan3::dna5;

    const seqan3::dna5_vector text1{"ACGACG"_dna5};
    const seqan3::dna5_vector text2{"ACGAACGC"_dna5};
    const seqan3::dna5_vector text3{"CGTCGT"_dna5};
    const seqan3::dna5_vector text4{"ATATAT"_dna5};
    const seqan3::dna5_vector empty_text{""_dna5};
    const size_t max_rank = seqan3::dna5::alphabet_size -1; // 'T'
};

// char - fm_index
template <typename sdsl_type>
struct fm_index_cursor_test<seqan3::fm_index_cursor<seqan3::fm_index<char,
                                                                    seqan3::text_layout::single,
                                                                    sdsl_type>>> : public ::testing::Test
{
    using alphabet_type = char;

    const std::string text1{"ACGACG"};
    const std::string text2{"ACGAACGC"};
    const std::string text3{"CGTCGT"};
    const std::string text4{"ATATAT"};
    const std::string empty_text{""};
    const size_t max_rank = CHAR_MAX;
};

// dna4
using it_t1 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(default_traits, fm_index_cursor_test, it_t1, );

using it_t2 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4,
                                                       seqan3::text_layout::single,
                                                       sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(byte_alphabet_traits, fm_index_cursor_test, it_t2, );

using it_t3 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_default_traits, fm_index_cursor_test, it_t3, );

using it_t4 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4,
                                                             seqan3::text_layout::single,
                                                             sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_byte_alphabet_traits, fm_index_cursor_test, it_t4, );

// dna5
using it_t5 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna5, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna5_default_traits, fm_index_cursor_test, it_t5, );

// char
using it_t6 = seqan3::fm_index_cursor<seqan3::fm_index<char, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(char_default_traits, fm_index_cursor_test, it_t6, );
