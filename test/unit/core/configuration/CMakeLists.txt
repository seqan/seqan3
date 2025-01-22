# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

seqan3_test (configuration_test.cpp)
seqan3_test (pipeable_config_element_test.cpp
             CYCLIC_DEPENDING_INCLUDES include-seqan3-core-configuration-configuration.hpp)
