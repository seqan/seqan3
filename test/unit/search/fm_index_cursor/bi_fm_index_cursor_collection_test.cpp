// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string_view>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/slice.hpp>

#include "bi_fm_index_cursor_collection_test_template.hpp"

// generic type
template <typename index_cursor_t>
struct bi_fm_index_cursor_collection_test : public ::testing::Test
{
    using index_type = typename index_cursor_t::index_type;
    using alphabet_type = typename index_type::alphabet_type;
    using text_type = std::vector<alphabet_type>;

    static constexpr auto convert = seqan3::views::char_to<alphabet_type> | seqan3::ranges::to<text_type>();

    text_type text{convert(std::string_view{"ACGGTAGGACGTAGC"})};
    text_type text1{convert(std::string_view{"AACGATCGGA"})};
    text_type text2{convert(std::string_view{"TGCTACGATCC"})};
    text_type text3 = seqan3::views::slice(text, 0, 11) | seqan3::ranges::to<text_type>(); // "ACGGTAGGACG"
    text_type text4 = seqan3::views::slice(text, 0, 14) | seqan3::ranges::to<text_type>(); // "ACGGTAGGACGTAG"

    std::vector<text_type> text_col1{text1, text1};
    std::vector<text_type> text_col2{text3, text2};
    std::vector<text_type> text_col3{text4, text2};
    std::vector<text_type> text_col4{text, text2};

    std::vector<text_type> rev_text1 =
        text_col1 | seqan3::views::deep{std::views::reverse} | seqan3::ranges::to<std::vector<text_type>>();
    std::vector<text_type> rev_text2 = text_col4 | seqan3::views::deep{std::views::reverse} | std::views::reverse
                                     | seqan3::ranges::to<std::vector<text_type>>();

    text_type pattern1{convert(std::string_view{"CAG"})};
    text_type pattern2{convert(std::string_view{"TT"})};
    text_type pattern3{convert(std::string_view{"GATGC"})};
    text_type pattern4{convert(std::string_view{"GATGG"})};
};

using it_t1 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, bi_fm_index_cursor_collection_test, it_t1, );

using it_t2 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, bi_fm_index_cursor_collection_test, it_t2, );
