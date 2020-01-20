// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <array>

#include <seqan3/core/algorithm/configuration_utility.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

enum class test_algo_id : uint8_t
{
    bar_id = 0,
    bax_id = 1,
    foo_id = 2,
    foobar_id = 3,
    SIZE = 4
};

namespace seqan3::detail
{

template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(::test_algo_id::SIZE)>,
                            static_cast<uint8_t>(::test_algo_id::SIZE)> compatibility_table<::test_algo_id>
{{
    {0, 1, 1, 1},
    {1, 0, 1, 0},
    {1, 1, 0, 1},
    {1, 0, 1, 0}
}};

}  // namespace seqan3::detail

struct bar : public seqan3::pipeable_config_element<bar, int>
{
    static constexpr test_algo_id id{test_algo_id::bar_id};
};

struct bax : public seqan3::pipeable_config_element<bax, float>
{
    static constexpr test_algo_id id{test_algo_id::bax_id};
};

struct foo : public seqan3::pipeable_config_element<foo, std::string>
{
    static constexpr test_algo_id id{test_algo_id::foo_id};
};

template <typename t = std::vector<int>>
struct foobar : public seqan3::pipeable_config_element<foobar<t>, t>
{
    static constexpr test_algo_id id{test_algo_id::foobar_id};
};
