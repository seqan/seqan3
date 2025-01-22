// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/aminoacid/aa27.hpp>

#include "fm_index_collection_test_template.hpp"
#include "fm_index_test_template.hpp"

using t1 = std::pair<seqan3::bi_fm_index<seqan3::aa27, seqan3::text_layout::single>, seqan3::aa27_vector>;
INSTANTIATE_TYPED_TEST_SUITE_P(aa27, fm_index_test, t1, );
using t2 =
    std::pair<seqan3::bi_fm_index<seqan3::aa27, seqan3::text_layout::collection>, std::vector<seqan3::aa27_vector>>;
INSTANTIATE_TYPED_TEST_SUITE_P(aa27_collection, fm_index_collection_test, t2, );
