// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

#include "../../core/algorithm/pipeable_config_element_test_template.hpp"

// ---------------------------------------------------------------------------------------------------------------------
// pipeable_config_element_test template
// ---------------------------------------------------------------------------------------------------------------------

using config_element_types = ::testing::Types<seqan3::detail::method_global_tag,
                                              seqan3::detail::method_local_tag>;

INSTANTIATE_TYPED_TEST_SUITE_P(method, pipeable_config_element_test, config_element_types, );
