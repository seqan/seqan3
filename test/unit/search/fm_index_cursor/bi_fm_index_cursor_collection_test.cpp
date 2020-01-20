// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "bi_fm_index_cursor_collection_test_template.hpp"

using it_t1 = bi_fm_index_cursor<bi_fm_index<dna4, text_layout::collection>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, bi_fm_index_cursor_collection_test, it_t1, );
