// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../edit_distance_unbanded_test_template.hpp"
#include "../fixture/global_edit_distance_unbanded.hpp"

using global_edit_distance_unbanded_types
    = ::testing::Types<
        global_fixture<&global::edit_distance::unbanded::dna4_01, uint8_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_01, uint16_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_01, uint32_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_01, uint64_t>,

        global_fixture<&global::edit_distance::unbanded::dna4_01T, uint8_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_01T, uint16_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_01T, uint32_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_01T, uint64_t>,

        global_fixture<&global::edit_distance::unbanded::dna4_02, uint8_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_02, uint16_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_02, uint32_t>,
        global_fixture<&global::edit_distance::unbanded::dna4_02, uint64_t>,

        global_fixture<&global::edit_distance::unbanded::aa27_01, uint8_t>,
        global_fixture<&global::edit_distance::unbanded::aa27_01, uint16_t>,
        global_fixture<&global::edit_distance::unbanded::aa27_01, uint32_t>,
        global_fixture<&global::edit_distance::unbanded::aa27_01, uint64_t>,

        global_fixture<&global::edit_distance::unbanded::aa27_01T, uint8_t>,
        global_fixture<&global::edit_distance::unbanded::aa27_01T, uint16_t>,
        global_fixture<&global::edit_distance::unbanded::aa27_01T, uint32_t>,
        global_fixture<&global::edit_distance::unbanded::aa27_01T, uint64_t>
    >;

INSTANTIATE_TYPED_TEST_CASE_P(global, edit_distance_unbanded, global_edit_distance_unbanded_types);
