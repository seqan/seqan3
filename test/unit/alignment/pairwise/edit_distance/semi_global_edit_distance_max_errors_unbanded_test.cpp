// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../edit_distance_unbanded_test_template.hpp"
#include "../fixture/semi_global_edit_distance_max_errors_unbanded.hpp"

using semi_global_edit_distance_max_errors_unbanded_types1
    = ::testing::Types<
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e5, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e5, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e5, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e5, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e2, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e2, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e2, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e2, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint64_t>
        >;

using semi_global_edit_distance_max_errors_unbanded_types2
    = ::testing::Types<
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e4, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e4, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e4, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e4, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e3, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e3, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e3, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e3, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e0, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e0, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e0, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s3u_15u_e0, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e0, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e0, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e0, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e0, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e5, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e5, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e5, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_s17u_1u_e5, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e0, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e0, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e0, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::dna4_03_e0, uint64_t>
    >;

using semi_global_edit_distance_max_errors_unbanded_types3
    = ::testing::Types<
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint64_t>,

        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint8_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint16_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint32_t>,
        semi_global_fixture<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint64_t>
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(semi_global_max_errors1, edit_distance_unbanded_test, semi_global_edit_distance_max_errors_unbanded_types1, );
INSTANTIATE_TYPED_TEST_SUITE_P(semi_global_max_errors2, edit_distance_unbanded_test, semi_global_edit_distance_max_errors_unbanded_types2, );
INSTANTIATE_TYPED_TEST_SUITE_P(semi_global_max_errors3, edit_distance_unbanded_test, semi_global_edit_distance_max_errors_unbanded_types3, );
