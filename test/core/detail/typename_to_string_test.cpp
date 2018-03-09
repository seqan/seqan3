// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#include <gtest/gtest.h>

#include <string>

#include <seqan3/core/detail/typename_to_string.hpp>

using namespace std::literals;
using namespace seqan3;

TEST(typename_to_string, builtin_types)
{
    EXPECT_EQ(detail::typename_to_string(static_cast<short>(0)),   "short"s);
    EXPECT_EQ(detail::typename_to_string(static_cast<int>(0)),     "int"s);
    EXPECT_EQ(detail::typename_to_string(static_cast<double>(0)),  "double"s);
    EXPECT_EQ(detail::typename_to_string(static_cast<uint8_t>(0)), "unsigned char"s);
    EXPECT_EQ(detail::typename_to_string("hello"),                 "char [6]"s);
    const char * val = "hello";
    EXPECT_EQ(detail::typename_to_string(val), "char const*"s);
}

template <typename type>
struct foo
{};

TEST(typename_to_string, flat_template_types)
{
    EXPECT_EQ(detail::typename_to_string(foo<short>{}),          "foo<short>"s);
    EXPECT_EQ(detail::typename_to_string(foo<int>{}),            "foo<int>"s);
    EXPECT_EQ(detail::typename_to_string(foo<double>{}),         "foo<double>"s);
    EXPECT_EQ(detail::typename_to_string(foo<uint8_t>{}),        "foo<unsigned char>"s);
    foo<decltype("hello")> f{};
    EXPECT_EQ(detail::typename_to_string(f), "foo<char const (&) [6]>"s);
}

template <typename type>
struct bar
{};

TEST(typename_to_string, nested_template_types)
{
    EXPECT_EQ(detail::typename_to_string(bar<foo<short>>{}),          "bar<foo<short> >"s);
    EXPECT_EQ(detail::typename_to_string(bar<foo<int>>{}),            "bar<foo<int> >"s);
    EXPECT_EQ(detail::typename_to_string(bar<foo<double>>{}),         "bar<foo<double> >"s);
    EXPECT_EQ(detail::typename_to_string(bar<foo<uint8_t>>{}),        "bar<foo<unsigned char> >"s);
    bar<foo<decltype("hello")>> f{};
    EXPECT_EQ(detail::typename_to_string(f), "bar<foo<char const (&) [6]> >"s);
}

template <typename param_t, template<typename> typename type>
struct foobar : bar<type<param_t>>
{};

TEST(typename_to_string, template_template_types)
{
    EXPECT_EQ((detail::typename_to_string(foobar<short, foo>{})),            "foobar<short, foo>"s);
    EXPECT_EQ((detail::typename_to_string(foobar<int, foo>{})),              "foobar<int, foo>"s);
    EXPECT_EQ((detail::typename_to_string(foobar<double, foo>{})),           "foobar<double, foo>"s);
    EXPECT_EQ((detail::typename_to_string(foobar<uint8_t, foo>{})),          "foobar<unsigned char, foo>"s);
    foobar<decltype("hello"), foo> f{};
    EXPECT_EQ((detail::typename_to_string(f)), "foobar<char const (&) [6], foo>"s);
}
