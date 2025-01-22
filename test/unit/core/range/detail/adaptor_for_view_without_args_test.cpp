// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <memory>

#include <seqan3/core/range/detail/adaptor_for_view_without_args.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

template <typename t>
struct dummy_view
{};

TEST(adaptor_combination, constexpr_combine)
{
    constexpr auto adaptor1 = seqan3::detail::adaptor_for_view_without_args<dummy_view>{};
    constexpr auto adaptor2 = seqan3::detail::adaptor_for_view_without_args<dummy_view>{};
    [[maybe_unused]] static constinit auto combined = adaptor1 | adaptor2;
}
