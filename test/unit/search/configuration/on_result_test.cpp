// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/search/configuration/on_result.hpp>

#include "../../core/configuration/pipeable_config_element_test_template.hpp"

// ---------------------------------------------------------------------------------------------------------------------
// test template : pipeable_config_element_test
// ---------------------------------------------------------------------------------------------------------------------

inline constexpr auto on_result_caller = [] (auto &&) {};

using test_types = ::testing::Types<seqan3::search_cfg::on_result<decltype(on_result_caller)>>;

INSTANTIATE_TYPED_TEST_SUITE_P(on_result_config, pipeable_config_element_test, test_types, );
