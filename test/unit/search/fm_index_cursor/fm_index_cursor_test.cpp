// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_cursor_test_template.hpp"

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
