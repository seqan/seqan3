// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest-spi.h> // provides test utility to test google test itself
#include <gtest/gtest.h>

#include <seqan3/test/expect_same_type.hpp>

TEST(expect_same_type, braces_with_many_commas)
{
    EXPECT_SAME_TYPE(decltype(std::tuple<int, float>{}), (std::tuple<int, float>));
    EXPECT_SAME_TYPE((decltype(std::tuple<int, float>{})), (std::tuple<int, float>));

    EXPECT_SAME_TYPE((std::tuple<int, float>), decltype(std::tuple<int, float>{}));
    EXPECT_SAME_TYPE((std::tuple<int, float>), (decltype(std::tuple<int, float>{})));
}

TEST(int, same_type_pass)
{
    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(0)>{}",
                                                             "std::type_identity< int>{}",
                                                             std::type_identity<decltype(0)>{},
                                                             std::type_identity<int>{});
    EXPECT_TRUE(expect_result);
    EXPECT_SAME_TYPE(decltype(0), int);
}

TEST(int, same_type_fail)
{
    char const * error_message = "Expected equality of these values:\n"
                                 "  decltype(0)\n"
                                 "    Which is: \"int\"\n"
                                 "  unsigned\n"
                                 "    Which is: \"unsigned int\"";

    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(0)>{}",
                                                             "std::type_identity< unsigned>{}",
                                                             std::type_identity<decltype(0)>{},
                                                             std::type_identity<unsigned>{});
    EXPECT_FALSE(expect_result);
    EXPECT_STREQ(error_message, expect_result.message());
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE(decltype(0), unsigned), error_message);

    // with parentheses has the same output
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE((decltype(0)), (unsigned)), error_message);
}

TEST(int_ref, same_type_pass)
{
    int a{};
    int & i = a;
    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(i)>{}",
                                                             "std::type_identity< int &>{}",
                                                             std::type_identity<decltype(i)>{},
                                                             std::type_identity<int &>{});
    EXPECT_TRUE(expect_result);
    EXPECT_SAME_TYPE(decltype(i), int &);
}

TEST(int_ref, same_type_fail)
{
    char const * error_message = "Expected equality of these values:\n"
                                 "  decltype(i)\n"
                                 "    Which is: \"int &\"\n"
                                 "  unsigned &\n"
                                 "    Which is: \"unsigned int &\"";

    int a{};
    int & i = a;
    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(i)>{}",
                                                             "std::type_identity< unsigned &>{}",
                                                             std::type_identity<decltype(i)>{},
                                                             std::type_identity<unsigned &>{});
    EXPECT_FALSE(expect_result);
    EXPECT_STREQ(error_message, expect_result.message());
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE(decltype(i), unsigned &), error_message);

    // with parentheses has the same output
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE((decltype(i)), (unsigned &)), error_message);
}

TEST(int_const_ref, same_type_pass)
{
    int const & i = 0;
    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(i)>{}",
                                                             "std::type_identity< int const &>{}",
                                                             std::type_identity<decltype(i)>{},
                                                             std::type_identity<int const &>{});
    EXPECT_TRUE(expect_result);
    EXPECT_SAME_TYPE(decltype(i), int const &);
}

TEST(int_const_ref, same_type_fail)
{
    char const * error_message = "Expected equality of these values:\n"
                                 "  decltype(i)\n"
                                 "    Which is: \"int const &\"\n"
                                 "  unsigned const &\n"
                                 "    Which is: \"unsigned int const &\"";

    int const & i = 0;
    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(i)>{}",
                                                             "std::type_identity< unsigned const &>{}",
                                                             std::type_identity<decltype(i)>{},
                                                             std::type_identity<unsigned const &>{});
    EXPECT_FALSE(expect_result);
    EXPECT_STREQ(error_message, expect_result.message());
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE(decltype(i), unsigned const &), error_message);

    // with parentheses has the same output
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE((decltype(i)), (unsigned const &)), error_message);
}

TEST(int_rvalue_ref, same_type_pass)
{
    int && i = 0;
    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(i)>{}",
                                                             "std::type_identity< int &&>{}",
                                                             std::type_identity<decltype(i)>{},
                                                             std::type_identity<int &&>{});
    EXPECT_TRUE(expect_result);
    EXPECT_SAME_TYPE(decltype(i), int &&);
}

TEST(int_rvalue_ref, same_type_fail)
{
    char const * error_message = "Expected equality of these values:\n"
                                 "  decltype(i)\n"
                                 "    Which is: \"int &&\"\n"
                                 "  unsigned &&\n"
                                 "    Which is: \"unsigned int &&\"";

    int && i = 0;
    auto && expect_result = seqan3::test::expect_same_type{}("std::type_identity< decltype(i)>{}",
                                                             "std::type_identity< unsigned &&>{}",
                                                             std::type_identity<decltype(i)>{},
                                                             std::type_identity<unsigned &&>{});
    EXPECT_FALSE(expect_result);
    EXPECT_STREQ(error_message, expect_result.message());
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE(decltype(i), unsigned &&), error_message);

    // with parentheses has the same output
    EXPECT_NONFATAL_FAILURE(EXPECT_SAME_TYPE((decltype(i)), (unsigned &&)), error_message);
}

TEST(tuple, same_type_pass)
{
    auto && expect_result =
        seqan3::test::expect_same_type{}("std::type_identity< decltype(std::tuple{0, .0f, .0, 0u})>{}",
                                         "std::type_identity< std::tuple<int, float, double, unsigned>>{}",
                                         std::type_identity<decltype(std::tuple{0, .0f, .0, 0u})>{},
                                         std::type_identity<std::tuple<int, float, double, unsigned>>{});

    EXPECT_TRUE(expect_result);
    EXPECT_SAME_TYPE(decltype(std::tuple{0, .0f, .0, 0u}), (std::tuple<int, float, double, unsigned>));
}

TEST(tuple, same_type_fail)
{
#ifdef _LIBCPP_VERSION
    char const * error_message = "Expected equality of these values:\n"
                                 "  decltype(std::tuple{0, .0f, .0, 0u})\n"
                                 "    Which is: \"std::__1::tuple<int, float, double, unsigned int>\"\n"
                                 "  std::tuple<int, float, unsigned, double>\n"
                                 "    Which is: \"std::__1::tuple<int, float, unsigned int, double>\"";
#else
    char const * error_message = "Expected equality of these values:\n"
                                 "  decltype(std::tuple{0, .0f, .0, 0u})\n"
                                 "    Which is: \"std::tuple<int, float, double, unsigned int>\"\n"
                                 "  std::tuple<int, float, unsigned, double>\n"
                                 "    Which is: \"std::tuple<int, float, unsigned int, double>\"";
#endif
    auto && expect_result =
        seqan3::test::expect_same_type{}("std::type_identity< decltype(std::tuple{0, .0f, .0, 0u})>{}",
                                         "std::type_identity< std::tuple<int, float, unsigned, double>>{}",
                                         std::type_identity<decltype(std::tuple{0, .0f, .0, 0u})>{},
                                         std::type_identity<std::tuple<int, float, unsigned, double>>{});

    EXPECT_FALSE(expect_result);
    EXPECT_STREQ(error_message, expect_result.message());
    EXPECT_NONFATAL_FAILURE(
        EXPECT_SAME_TYPE(decltype(std::tuple{0, .0f, .0, 0u}), (std::tuple<int, float, unsigned, double>)),
        error_message);
}
