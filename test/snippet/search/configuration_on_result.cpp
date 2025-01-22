// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/configuration/on_result.hpp>

int main()
{
    seqan3::search_cfg::on_result cfg{[](auto && result)
                                      {
                                          seqan3::debug_stream << result << '\n';
                                      }};
    return 0;
}
