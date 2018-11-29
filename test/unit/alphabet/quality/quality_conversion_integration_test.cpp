// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <numeric>
#include <sstream>
#include <tuple>
#include <utility>

#include <gtest/gtest.h>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/quality/all.hpp>

using namespace seqan3;

template <typename T>
class quality_conversion : public ::testing::Test
{};

// add all alphabets from the quality sub module here
using quality_conversion_types     = ::testing::Types<phred42, phred63, phred68legacy>;

TYPED_TEST_CASE(quality_conversion, quality_conversion_types);

template<typename Function, typename... Args>
void va_for_each(Function&& fn, Args&&... args)
{
    using aux = int[];
    static_cast<void>(aux
        {
            0, (static_cast<void>(fn(std::forward<Args>(args))), 0)...
        });
}

TYPED_TEST(quality_conversion, implicit_conversion)
{
    auto convert = [](auto other) {
        if constexpr (!std::is_same_v<TypeParam, decltype(other)>)
        {
            TypeParam p{0};
            p.assign_phred(0);
            [[maybe_unused]] decltype(other) p_other = static_cast<decltype(other)>(p);
            EXPECT_EQ(p_other.to_phred(), typename decltype(other)::phred_type(0));
        }
    };
    phred42 p42;
    phred63 p63;
    phred68legacy p68;
    va_for_each(convert, p42, p63, p68);
}
