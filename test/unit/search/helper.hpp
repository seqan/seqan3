// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <iterator>
#include <vector>

#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#include <seqan3/search/fm_index/fm_index_cursor.hpp>
#include <seqan3/utility/range/to.hpp>

namespace seqan3
{
template <typename char_t, typename index_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, seqan3::fm_index_cursor<index_t> const &)
{
    return s << ("fm_index_cursor");
}

template <typename char_t, typename index_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s,
                                              seqan3::bi_fm_index_cursor<index_t> const &)
{
    return s << ("bi_fm_index_cursor");
}

template <typename result_range_t>
std::vector<std::ranges::range_value_t<result_range_t>> uniquify(result_range_t && result_range)
{
    auto unique_res = result_range | ranges::to<std::vector>();
    std::sort(unique_res.begin(), unique_res.end());
    unique_res.erase(std::unique(unique_res.begin(), unique_res.end()), unique_res.end());
    return unique_res;
}

} // namespace seqan3
