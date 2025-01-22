// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "../edit_distance_unbanded_test_template.hpp"
#include "../fixture/global_edit_distance_unbanded.hpp"

using global_edit_distance_unbanded_types =
    ::testing::Types<global_fixture<&global::edit_distance::unbanded::dna4_01, uint8_t>,
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

                     global_fixture<&global::edit_distance::unbanded::dna4_02_s10u_15u, uint8_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s10u_15u, uint16_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s10u_15u, uint32_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s10u_15u, uint64_t>,

                     global_fixture<&global::edit_distance::unbanded::dna4_02_s3u_15u, uint8_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s3u_15u, uint16_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s3u_15u, uint32_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s3u_15u, uint64_t>,

                     global_fixture<&global::edit_distance::unbanded::dna4_02_s1u_15u, uint8_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s1u_15u, uint16_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s1u_15u, uint32_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02_s1u_15u, uint64_t>,

                     global_fixture<&global::edit_distance::unbanded::dna4_02T_s15u_1u, uint8_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02T_s15u_1u, uint16_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02T_s15u_1u, uint32_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_02T_s15u_1u, uint64_t>,

                     global_fixture<&global::edit_distance::unbanded::dna4_03, uint8_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_03, uint16_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_03, uint32_t>,
                     global_fixture<&global::edit_distance::unbanded::dna4_03, uint64_t>,

                     global_fixture<&global::edit_distance::unbanded::aa27_01, uint8_t>,
                     global_fixture<&global::edit_distance::unbanded::aa27_01, uint16_t>,
                     global_fixture<&global::edit_distance::unbanded::aa27_01, uint32_t>,
                     global_fixture<&global::edit_distance::unbanded::aa27_01, uint64_t>,

                     global_fixture<&global::edit_distance::unbanded::aa27_01T, uint8_t>,
                     global_fixture<&global::edit_distance::unbanded::aa27_01T, uint16_t>,
                     global_fixture<&global::edit_distance::unbanded::aa27_01T, uint32_t>,
                     global_fixture<&global::edit_distance::unbanded::aa27_01T, uint64_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(global, edit_distance_unbanded_test, global_edit_distance_unbanded_types, );
