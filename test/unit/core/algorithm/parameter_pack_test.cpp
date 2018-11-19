// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

TEST(parameter_pack, for_each_value)
{
    int i = 0;
    auto fn = [&i](int arg)
    {
        EXPECT_EQ(i, arg);
        ++i;
    };

    detail::for_each_value(fn);
    EXPECT_EQ(i, 0);
    detail::for_each_value(fn, 0);
    EXPECT_EQ(i, 1);
    detail::for_each_value(fn, 1, 2);
    EXPECT_EQ(i, 3);
    detail::for_each_value(fn, 3, 4, 5);
    EXPECT_EQ(i, 6);
}

TEST(parameter_pack, for_each_value2)
{
    std::stringstream s{};
    debug_stream_type stream{s};

    auto fn = [&stream](auto && arg)
    {
        stream << arg << ";";
    };

    detail::for_each_value(fn);
    detail::for_each_value(fn, 0);
    detail::for_each_value(fn, 1.0, '2');
    detail::for_each_value(fn, "3;4", -5, 'C'_dna4);

    EXPECT_EQ(s.str(), "0;1;2;3;4;-5;C;");
}

TEST(parameter_pack, for_each_type)
{
    std::stringstream s{};
    debug_stream_type stream{s};

    auto fn = [&stream](auto id)
    {
        using type = typename decltype(id)::type;

        // id is of type std::type_identity<type>
        using id_t = decltype(id);
        static_assert(std::is_same_v<id_t, std::type_identity<type>>);

        if constexpr (std::is_same_v<type, bool>)
            stream << type{0} << ";";
        if constexpr (std::is_same_v<type, uint8_t>)
            stream << type{1} << ";";
        if constexpr (std::is_same_v<type, int8_t>)
            stream << type{-1} << ";";
        if constexpr (std::is_same_v<type, uint16_t>)
            stream << type{2} << ";";
        if constexpr (std::is_same_v<type, int16_t>)
            stream << type{-2} << ";";
        if constexpr (std::is_same_v<type, uint32_t>)
            stream << type{3} << ";";
        if constexpr (std::is_same_v<type, int32_t>)
            stream << type{-3} << ";";
        if constexpr (std::is_same_v<type, uint64_t>)
            stream << type{4} << ";";
        if constexpr (std::is_same_v<type, int64_t>)
            stream << type{-4} << ";";
    };

    detail::for_each_type<bool, uint8_t, int8_t, uint16_t, int16_t,
                          uint32_t, int32_t, uint64_t, int64_t>(fn);

    EXPECT_EQ(s.str(), "0;1;-1;2;-2;3;-3;4;-4;");
}

TEST(type_list, for_each_type)
{
    std::stringstream s{};
    debug_stream_type stream{s};

    auto fn = [&stream](auto id)
    {
        using type = typename decltype(id)::type;

        // id is of type std::type_identity<type>
        using id_t = decltype(id);
        static_assert(std::is_same_v<id_t, std::type_identity<type>>);

        if constexpr (std::is_same_v<type, bool>)
            stream << type{0} << ";";
        if constexpr (std::is_same_v<type, uint8_t>)
            stream << type{1} << ";";
        if constexpr (std::is_same_v<type, int8_t>)
            stream << type{-1} << ";";
        if constexpr (std::is_same_v<type, uint16_t>)
            stream << type{2} << ";";
        if constexpr (std::is_same_v<type, int16_t>)
            stream << type{-2} << ";";
        if constexpr (std::is_same_v<type, uint32_t>)
            stream << type{3} << ";";
        if constexpr (std::is_same_v<type, int32_t>)
            stream << type{-3} << ";";
        if constexpr (std::is_same_v<type, uint64_t>)
            stream << type{4} << ";";
        if constexpr (std::is_same_v<type, int64_t>)
            stream << type{-4} << ";";
    };

    using types = type_list<bool, uint8_t, int8_t, uint16_t, int16_t,
                            uint32_t, int32_t, uint64_t, int64_t>;
    detail::for_each_type(fn, types{});

    EXPECT_EQ(s.str(), "0;1;-1;2;-2;3;-3;4;-4;");
}

TEST(tuple, for_each_type)
{
    std::stringstream s{};
    debug_stream_type stream{s};

    auto fn = [&stream](auto id)
    {
        using type = typename decltype(id)::type;

        // id is of type std::type_identity<type>
        using id_t = decltype(id);
        static_assert(std::is_same_v<id_t, std::type_identity<type>>);

        if constexpr (std::is_same_v<type, bool>)
            stream << type{0} << ";";
        if constexpr (std::is_same_v<type, uint8_t>)
            stream << type{1} << ";";
        if constexpr (std::is_same_v<type, int8_t>)
            stream << type{-1} << ";";
        if constexpr (std::is_same_v<type, uint16_t>)
            stream << type{2} << ";";
        if constexpr (std::is_same_v<type, int16_t>)
            stream << type{-2} << ";";
        if constexpr (std::is_same_v<type, uint32_t>)
            stream << type{3} << ";";
        if constexpr (std::is_same_v<type, int32_t>)
            stream << type{-3} << ";";
        if constexpr (std::is_same_v<type, uint64_t>)
            stream << type{4} << ";";
        if constexpr (std::is_same_v<type, int64_t>)
            stream << type{-4} << ";";
    };

    using tuple = std::tuple<bool, uint8_t, int8_t, uint16_t, int16_t,
                             uint32_t, int32_t, uint64_t, int64_t>;
    detail::for_each_type(fn, tuple{});

    EXPECT_EQ(s.str(), "0;1;-1;2;-2;3;-3;4;-4;");
}
