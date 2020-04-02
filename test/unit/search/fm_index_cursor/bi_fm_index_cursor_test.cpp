// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <string>
#include <string_view>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/ranges>

#include "bi_fm_index_cursor_test_template.hpp"

// generic type
template <typename index_cursor_t>
struct bi_fm_index_cursor_test : public ::testing::Test
{
    using index_type = typename index_cursor_t::index_type;
    using alphabet_type = typename index_type::alphabet_type;
    using text_type = std::vector<alphabet_type>;

    static constexpr auto convert = seqan3::views::char_to<alphabet_type> | seqan3::views::to<text_type>;

    text_type text{convert(std::string_view{"ACGGTAGGACGTAGC"})};
    text_type text1{convert(std::string_view{"AACGATCGGA"})};

    text_type rev_text1 = text | std::views::reverse | seqan3::views::to<text_type>;
    text_type rev_text2 = text1 | std::views::reverse | seqan3::views::to<text_type>;

    text_type pattern1{convert(std::string_view{"CAG"})};
    text_type pattern2{convert(std::string_view{"TT"})};
    text_type pattern3{convert(std::string_view{"GATGC"})};
    text_type pattern4{convert(std::string_view{"GATGG"})};
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
