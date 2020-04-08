// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/std/concepts>

TEST(traits, alignment_function_traits)
{
    using result_callback_t = std::function<void(int)>;  // The callback type.
    using algorithm_t = std::function<void(std::string, result_callback_t)>; // The algorithm type.

    // Get the traits class and evaluate the member types.
    using function_traits_t =  seqan3::detail::alignment_function_traits<algorithm_t>;

    using input_t = typename function_traits_t::sequence_input_type;
    using callback_t = typename function_traits_t::callback_type;
    using alignment_result_t = typename function_traits_t::alignment_result_type;

    EXPECT_TRUE((std::same_as<input_t, std::string>));
    EXPECT_TRUE((std::same_as<callback_t, result_callback_t>));
    EXPECT_TRUE((std::same_as<alignment_result_t, int>));
}
