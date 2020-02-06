// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/builtin_simd.hpp>
#include <seqan3/core/simd/detail/default_simd_backend.hpp>

#include <iostream>
#include <type_traits>

TEST(default_simd_backend, test)
{
    EXPECT_TRUE((std::is_same_v<seqan3::detail::default_simd_backend<int16_t, 8>,
                                seqan3::detail::builtin_simd<int16_t, 8>>));
}
