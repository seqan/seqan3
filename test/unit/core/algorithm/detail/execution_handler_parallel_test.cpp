// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>

#include "execution_handler_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(execution_handler_parallel,
                               execution_handler,
                               seqan3::detail::execution_handler_parallel, );
