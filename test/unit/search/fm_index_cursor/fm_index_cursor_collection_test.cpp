// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_cursor_collection_test_template.hpp"

using it_t1 = fm_index_cursor<fm_index<dna4, text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(default_traits, fm_index_cursor_collection_test, it_t1, );

using it_t2 = fm_index_cursor<fm_index<dna4, text_layout::collection, sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(byte_alphabet_traits, fm_index_cursor_collection_test, it_t2, );

using it_t3 = bi_fm_index_cursor<bi_fm_index<dna4, text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_default_traits, fm_index_cursor_collection_test, it_t3, );

using it_t4 = bi_fm_index_cursor<bi_fm_index<dna4, text_layout::collection, sdsl_byte_index_type>>;
INSTANTIATE_TYPED_TEST_SUITE_P(bi_byte_alphabet_traits, fm_index_cursor_collection_test, it_t4, );
