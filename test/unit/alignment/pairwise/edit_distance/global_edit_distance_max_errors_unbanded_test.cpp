// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "../edit_distance_unbanded_test_template.hpp"
#include "../fixture/global_edit_distance_max_errors_unbanded.hpp"

using global_edit_distance_max_errors_unbanded_types1 =
    ::testing::Types<global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e8, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e8, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e8, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e8, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e7, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e7, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e7, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e7, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e5, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e5, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e5, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01_e5, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint64_t>>;

using global_edit_distance_max_errors_unbanded_types2 =
    ::testing::Types<global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e8, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e8, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e8, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e8, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e4, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e4, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e4, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_e4, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e7, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e7, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e7, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s10u_15u_e7, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e255, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e5, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e5, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e5, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02_s1u_15u_e5, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e255, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e5, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e5, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e5, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_02T_s15u_1u_e5, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_03_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_03_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_03_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::dna4_03_e255, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint64_t>,

                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint8_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint16_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint32_t>,
                     global_fixture<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint64_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(global1, edit_distance_unbanded_test, global_edit_distance_max_errors_unbanded_types1, );
INSTANTIATE_TYPED_TEST_SUITE_P(global2, edit_distance_unbanded_test, global_edit_distance_max_errors_unbanded_types2, );
