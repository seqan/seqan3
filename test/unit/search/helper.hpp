// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <vector>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/iterator>

namespace seqan3
{

template <typename search_result_range_t>
std::vector<std::ranges::range_value_t<search_result_range_t>> uniquify(search_result_range_t && search_result_range)
{
    auto unique_res = search_result_range | views::to<std::vector>;
    std::sort(unique_res.begin(), unique_res.end());
    unique_res.erase(std::unique(unique_res.begin(), unique_res.end()), unique_res.end());
    return unique_res;
}

} // namespace std
