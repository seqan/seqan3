// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/configuration/configuration.hpp>

enum struct my_id : int
{
    bar_id,
    foo_id
};

namespace seqan3::detail
{
template <>
inline constexpr std::array<std::array<int, 2>, 2> compatibility_table<my_id>{{{0, 1}, {1, 0}}};
} // namespace seqan3::detail
