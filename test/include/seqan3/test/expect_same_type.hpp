// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides EXPECT_SAME_TYPE.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/utility/detail/type_name_as_string.hpp>

namespace seqan3::test
{
// https://stackoverflow.com/a/62984543
#define EXPECT_SAME_TYPE_DEPAREN(X) EXPECT_SAME_TYPE_ESC(EXPECT_SAME_TYPE_ISH X)
#define EXPECT_SAME_TYPE_ISH(...) EXPECT_SAME_TYPE_ISH __VA_ARGS__
#define EXPECT_SAME_TYPE_ESC(...) EXPECT_SAME_TYPE_ESC_(__VA_ARGS__)
#define EXPECT_SAME_TYPE_ESC_(...) EXPECT_SAME_TYPE_VAN##__VA_ARGS__
#define EXPECT_SAME_TYPE_VANEXPECT_SAME_TYPE_ISH

#define EXPECT_SAME_TYPE(val1, val2)                                                                                   \
    EXPECT_PRED_FORMAT2(::seqan3::test::expect_same_type{},                                                            \
                        (std::type_identity<EXPECT_SAME_TYPE_DEPAREN(val1)>{}),                                        \
                        (std::type_identity<EXPECT_SAME_TYPE_DEPAREN(val2)>{}));

struct expect_same_type
{
    template <typename lhs_t, typename rhs_t>
    ::testing::AssertionResult operator()(std::string lhs_expression,
                                          std::string rhs_expression,
                                          std::type_identity<lhs_t>,
                                          std::type_identity<rhs_t>)
    {
        auto remove_wrap_type_identity = [](std::string str)
        {
            // EXPECT_SAME_TYPE_DEPAREN adds a space after the prefix
            std::string prefix = "std::type_identity< ";
            std::string suffix = ">{}";

            size_t str_start = str.find(prefix) + prefix.size();
            size_t str_end = str.rfind(suffix);
            assert(str_end >= str_start);

            return str.substr(str_start, str_end - str_start);
        };

        if (std::is_same_v<lhs_t, rhs_t>)
            return ::testing::AssertionSuccess();

        return ::testing::internal::CmpHelperEQFailure(remove_wrap_type_identity(lhs_expression).c_str(),
                                                       remove_wrap_type_identity(rhs_expression).c_str(),
                                                       seqan3::detail::type_name_as_string<lhs_t>,
                                                       seqan3::detail::type_name_as_string<rhs_t>);
    }
};

} // namespace seqan3::test
