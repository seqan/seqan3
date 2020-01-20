// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "fm_index_collection_test_template.hpp"
#include "fm_index_test_template.hpp"

using t1 = std::pair<bi_fm_index<dna4, text_layout::single>, std::vector<dna4>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, fm_index_test, t1, );
using t2 = std::pair<bi_fm_index<dna4, text_layout::collection>, std::vector<std::vector<dna4>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(dna4_collection, fm_index_collection_test, t2, );
