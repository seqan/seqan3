// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <climits>
#include <string>
#include <string_view>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/char_to.hpp>

#include "fm_index_cursor_collection_test_template.hpp"

// generic type
template <typename index_cursor_t>
struct fm_index_cursor_collection_test : public ::testing::Test
{
    using index_type = typename index_cursor_t::index_type;
    using alphabet_type = typename index_type::alphabet_type;
    using text_type = std::vector<alphabet_type>;

    static constexpr auto convert = seqan3::views::char_to<alphabet_type> | seqan3::views::to<text_type>;

    text_type text1{convert(std::string_view{"ACGACG"})};
    text_type text2{convert(std::string_view{"ACGAACGC"})};
    text_type text3{convert(std::string_view{"CGTCGT"})};
    text_type text4{convert(std::string_view{"ATATAT"})};
    text_type text5{convert(std::string_view{"TGCGATCGA"})};
    text_type text6{convert(std::string_view{"TACGATCGA"})};
    text_type text7{convert(std::string_view{"ACGTACGT"})};
    text_type text8{convert(std::string_view{"TGCGATACGA"})};

    text_type empty_text{};

    std::vector<text_type> text_col1{text1, text1};
    std::vector<text_type> text_col2{text1, text5};
    std::vector<text_type> text_col3{text1, empty_text, empty_text, text5};
    std::vector<text_type> text_col4{text2, text6};
    std::vector<text_type> text_col5{text2, text5};
    std::vector<text_type> text_col6{text3, text3};
    std::vector<text_type> text_col7{text4, text4};
    std::vector<text_type> text_col8{text7, text8};
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
