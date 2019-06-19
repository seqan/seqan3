// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_cursor_test_template.hpp"

using it_t1 = fm_index_cursor<fm_index<text_layout::single>>;
INSTANTIATE_TYPED_TEST_CASE_P(default_traits, fm_index_cursor_test, it_t1);

using it_t2 = fm_index_cursor<fm_index<text_layout::single, sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_CASE_P(byte_alphabet_traits, fm_index_cursor_test, it_t2);

using it_t3 = bi_fm_index_cursor<bi_fm_index<text_layout::single>>;
INSTANTIATE_TYPED_TEST_CASE_P(bi_default_traits, fm_index_cursor_test, it_t3);

using it_t4 = bi_fm_index_cursor<bi_fm_index<text_layout::single, sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_CASE_P(bi_byte_alphabet_traits, fm_index_cursor_test, it_t4);
