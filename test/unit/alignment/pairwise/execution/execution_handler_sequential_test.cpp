// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/execution/execution_handler_sequential.hpp>

#include "execution_handler_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(execution_handler_sequential,
                               execution_handler,
                               seqan3::detail::execution_handler_sequential, );
