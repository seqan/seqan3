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

namespace seqan3
{

template <typename T>
std::vector<T> uniquify(std::vector<T> v)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

template <typename T>
std::vector<std::vector<T>> uniquify(std::vector<std::vector<T>> v)
{
    std::for_each(v.begin(), v.end(), [] (auto & hits) { uniquify(hits); } );
    return v;
}

template <typename T1, typename T2>
std::vector<std::pair<T1, T2>> uniquify(std::vector<std::pair<T1, T2>> v)
{
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

template <typename T1, typename T2>
std::vector<std::vector<std::pair<T1, T2>>> uniquify(std::vector<std::vector<std::pair<T1, T2>>> v)
{
    std::for_each(v.begin(), v.end(), [] (auto & hits) { uniquify(hits); } );
    return v;
}

void random_text(std::vector<dna4> & text, uint64_t const length)
{
    uint8_t alphabet_size{4};

    text.resize(length);
    for (uint64_t i = 0; i < length; ++i)
        assign_rank_to(std::rand() % alphabet_size, text[i]);
}

} // namespace std
