// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <string>
#include <string_view>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/fm_index_cursor.hpp>

#include "fm_index_cursor_test_template.hpp"

// generic type
template <typename index_cursor_t>
struct fm_index_cursor_test : public ::testing::Test
{
    using index_type = typename index_cursor_t::index_type;
    using alphabet_type = typename index_type::alphabet_type;
    using text_type = std::vector<alphabet_type>;

    static constexpr bool is_bi_fm_index =
        seqan3::detail::template_specialisation_of<typename index_type::cursor_type, seqan3::bi_fm_index_cursor>;
    static constexpr auto convert = seqan3::views::char_to<alphabet_type> | seqan3::ranges::to<text_type>();

    text_type text1{convert(std::string_view{"ACGACG"})};
    text_type text2{convert(std::string_view{"ACGAACGC"})};
    text_type text3{convert(std::string_view{"CGTCGT"})};
    text_type text4{convert(std::string_view{"ATATAT"})};
    text_type empty_text{};
};

// dna4
using it_t1 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(default_traits, fm_index_cursor_test, it_t1, );

using it_t2 =
    seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single, sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(byte_alphabet_traits, fm_index_cursor_test, it_t2, );

using it_t3 = seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_default_traits, fm_index_cursor_test, it_t3, );

using it_t4 =
    seqan3::bi_fm_index_cursor<seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single, sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_byte_alphabet_traits, fm_index_cursor_test, it_t4, );

// dna5
using it_t5 = seqan3::fm_index_cursor<seqan3::fm_index<seqan3::dna5, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna5_default_traits, fm_index_cursor_test, it_t5, );

// char
using it_t6 = seqan3::fm_index_cursor<seqan3::fm_index<char, seqan3::text_layout::single>>;
INSTANTIATE_TYPED_TEST_SUITE_P(char_default_traits, fm_index_cursor_test, it_t6, );
