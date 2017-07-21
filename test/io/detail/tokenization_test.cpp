// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <gtest/gtest.h>

#include <string>
#include <vector>
#include <iterator>
#include <type_traits>

#include <range/v3/utility/associated_types.hpp>
#include <range/v3/utility/iterator.hpp>
#include <range/v3/istream_range.hpp>
#include <range/v3/view/slice.hpp>

#include <seqan3/io/detail/tokenization.hpp>

#include "test_small_stream_buffer.hpp"

using namespace seqan3;

template <typename T>
class tokenization_test : public ::testing::Test
{};

// add all alphabets from the aminoacid sub module here
using source_target_types = ::testing::Types<std::pair<std::string, std::string>,
                                             std::pair<std::string, std::ostringstream>,
                                             std::pair<std::istringstream, std::string>,
                                             std::pair<std::istringstream, std::ostringstream>>;

TYPED_TEST_CASE(tokenization_test, source_target_types);
// What we actually want is a direction_iterator over stl containers and
// And we want to make them applicable to streams

// For example our own stream iterator can do this, if we pass our own stream_adapter.
// We want them chunkable if the container allows for it.

template <typename t>
auto test_tokenization_get_iterator(t && obj)
{
    if constexpr (std::is_same_v<std::decay_t<t>, std::string>)
        return std::begin(obj);
    else if constexpr (std::is_same_v<std::decay_t<t>, std::istringstream>)
        return std::istreambuf_iterator<ranges::value_type_t<std::decay_t<t>>>(obj);
    else
        return ranges::ostream_iterator<ranges::value_type_t<std::decay_t<t>>>(obj);
}

template <typename t>
auto test_tokenization_get_end([[maybe_unused]] t && obj)
{
    if constexpr (std::is_same_v<std::decay_t<t>, std::string>)
        return std::end(obj);
    else if constexpr (std::is_same_v<std::decay_t<t>, std::istringstream>)
        return std::istreambuf_iterator<ranges::value_type_t<std::decay_t<t>>>{};
}

template <typename t>
auto test_tokenization_get_data(t && obj)
{
    if constexpr (std::is_same_v<std::decay_t<t>, std::string>)
        return obj;
    else
        return obj.str();
}

TYPED_TEST(tokenization_test, write)
{
    using namespace seqan3::detail;

    using source_type = typename std::tuple_element<0, TypeParam>::type;
    using target_type = typename std::tuple_element<1, TypeParam>::type;

    {  // standard input + standard output
        source_type in{"hello_world"};
        target_type out{};
        write(test_tokenization_get_iterator(in), 11, out);
        EXPECT_EQ(test_tokenization_get_data(out), test_tokenization_get_data(in));
    }

    {  // standard input + preferred output iterator
        source_type in{"hello_world"};
        target_type out{};
        write(test_tokenization_get_iterator(in), 11, make_preferred_output_iterator(out));
        EXPECT_EQ(test_tokenization_get_data(out), test_tokenization_get_data(in));
    }

    {  // preferred input iterator + standard output
        source_type in{"hello_world"};
        target_type out{};
        auto rng = make_preferred_input_iterator_range(in);
        write(std::get<0>(rng), 11, out);
        EXPECT_EQ(test_tokenization_get_data(out), test_tokenization_get_data(in));
    }

    {  // preferred input iterator + preferred output iterator
        source_type in{"hello_world"};
        target_type out{};
        auto rng = make_preferred_input_iterator_range(in);
        write(std::get<0>(rng), 11, make_preferred_output_iterator(out));
        EXPECT_EQ(test_tokenization_get_data(out), test_tokenization_get_data(in));
    }

    if constexpr (std::is_same_v<source_type, std::string>)
    {  // container short cut.
        source_type in{"hello_world"};
        target_type out{};
        write(in, out);
        EXPECT_EQ(test_tokenization_get_data(out), test_tokenization_get_data(in));
    }
}

TYPED_TEST(tokenization_test, read_impl)
{
    using namespace seqan3::detail;

    using source_type = typename std::tuple_element<0, TypeParam>::type;
    using target_type = typename std::tuple_element<1, TypeParam>::type;

    {  // istream interface.
        source_type in{"hello_world"};
        target_type out;

        auto it = test_tokenization_get_iterator(in);
        auto o_iter = make_preferred_output_iterator(out);
        read_impl(it, test_tokenization_get_end(in), o_iter, equals_char<'_'>{}, equals_char<'o'>{});
        EXPECT_EQ(test_tokenization_get_data(out), "hell");
        read_impl(it, test_tokenization_get_end(in), o_iter, equals_char<'\n'>{}, equals_char<'l'>{});
        EXPECT_EQ(test_tokenization_get_data(out), "hell_word");
    }
}

template <typename in_t, typename sentinel_t, typename out_t, typename target_t>
void impl_tokenization_test_read_impl_chunked(in_t & in_cur, sentinel_t && in_end, out_t & out_it, target_t & target)
{
    read_impl(in_cur, in_end, out_it, seqan3::detail::equals_char<'_'>{}, seqan3::detail::equals_char<'o'>{});
    EXPECT_EQ(static_cast<std::string>(target | ranges::view::slice(0, 4)), "hell");
    read_impl(in_cur, in_end, out_it, seqan3::detail::equals_char<'\n'>{}, seqan3::detail::equals_char<'l'>{});
    EXPECT_EQ(static_cast<std::string>(target | ranges::view::slice(0, 9)), "hell_word");
}

TYPED_TEST(tokenization_test, read_impl_chunked)
{
    using namespace seqan3::detail;

    using source_type = typename std::tuple_element<0, TypeParam>::type;
    using target_type = std::conditional_t<std::is_same_v<typename std::tuple_element<1, TypeParam>::type, std::ostringstream>,
                                           std::ostream,
                                           typename std::tuple_element<1, TypeParam>::type>;

    if constexpr (std::is_same_v<target_type, std::ostream>)
    {  // istream interface.
        source_type in{"hello_world"};

        std::string storage;
        storage.resize(11);
        io_test_small_stream_buffer buf(storage.data(), storage.data() + storage.size());
        target_type out{&buf};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_impl_chunked(r_beg, r_end, o_iter, storage);
    }
    else
    {  // istream interface.
        source_type in{"hello_world"};
        target_type out{};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_impl_chunked(r_beg, r_end, o_iter, out);
    }
}

TYPED_TEST(tokenization_test, ignore_impl)
{
    using namespace seqan3::detail;

    using source_type = typename std::tuple_element<0, TypeParam>::type;

    {  // istream interface.
        source_type in{"hello_world"};

        auto it = test_tokenization_get_iterator(in);
        ignore_impl(it, test_tokenization_get_end(in), equals_char<'_'>{});
        EXPECT_EQ(*it, '_');
        ignore_impl(it, test_tokenization_get_end(in), equals_char<'d'>{});
        EXPECT_EQ(*it, 'd');
    }
}

TYPED_TEST(tokenization_test, ignore_impl_chunked)
{
    using namespace seqan3::detail;

    using source_type = typename std::tuple_element<0, TypeParam>::type;

    {  // istream interface.
        source_type in{"hello_world"};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        ignore_impl(r_beg, r_end, equals_char<'_'>{});
        EXPECT_EQ(*r_beg, '_');
        ignore_impl(r_beg, r_end, equals_char<'d'>{});
        EXPECT_EQ(*r_beg, 'd');
    }
}

template <typename in_t, typename sentinel_t, typename out_t, typename target_t>
void impl_tokenization_test_read_until(in_t & in_cur, sentinel_t && in_end, out_t & out_it, target_t & target)
{
    read_until(in_cur, in_end, out_it, seqan3::detail::equals_char<'_'>{}, seqan3::detail::equals_char<'o'>{});
    auto res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 4)), "hell");
    EXPECT_EQ(*in_cur, '_');
    read_until(in_cur, in_end, out_it, seqan3::detail::is_newline{});
    res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 10)), "hell_world");
    EXPECT_EQ(*in_cur, '\n');
}

TYPED_TEST(tokenization_test, read_until)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;
    using target_type   = typename std::tuple_element<1, TypeParam>::type;
    using target_type_2 = std::conditional_t<std::is_same_v<target_type, std::ostringstream>,
                                             std::ostream, target_type>;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        target_type out{};
        auto it = test_tokenization_get_iterator(in);
        auto o_iter = make_preferred_output_iterator(out);

        impl_tokenization_test_read_until(it, test_tokenization_get_end(in), o_iter, out);
    }

    if constexpr (std::is_same_v<target_type_2, std::ostream>)
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        std::string storage;
        storage.resize(11);
        io_test_small_stream_buffer buf(storage.data(), storage.data() + storage.size());
        target_type_2 out{&buf};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_until(r_beg, r_end, o_iter, storage);
    }
    else
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        target_type_2 out{};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_until(r_beg, r_end, o_iter, out);
    }
}

template <typename in_t, typename sentinel_t, typename out_t, typename target_t>
void impl_tokenization_test_read_n(in_t & in_cur, sentinel_t && in_end, out_t & out_it, target_t & target)
{
    read_n(in_cur, in_end, out_it, 5, seqan3::detail::equals_char<'o'>{});
    auto res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 4)), "hell");
    EXPECT_EQ(*in_cur, '_');

    read_n(in_cur, in_end, out_it, 6);
    res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 10)), "hell_world");
    EXPECT_EQ(*in_cur, '\n');
}

TYPED_TEST(tokenization_test, read_n)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;
    using target_type   = typename std::tuple_element<1, TypeParam>::type;
    using target_type_2 = std::conditional_t<std::is_same_v<target_type, std::ostringstream>,
                                             std::ostream, target_type>;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        target_type out{};
        auto it = test_tokenization_get_iterator(in);
        auto o_iter = make_preferred_output_iterator(out);

        impl_tokenization_test_read_n(it, test_tokenization_get_end(in), o_iter, out);
    }

    if constexpr (std::is_same_v<target_type_2, std::ostream>)
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        std::string storage;
        storage.resize(11);
        io_test_small_stream_buffer buf(storage.data(), storage.data() + storage.size());
        target_type_2 out{&buf};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_n(r_beg, r_end, o_iter, storage);
    }
    else
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        target_type_2 out{};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_n(r_beg, r_end, o_iter, out);
    }
}

template <typename in_t, typename sentinel_t, typename out_t, typename target_t>
void impl_tokenization_test_read_one(in_t & in_cur, sentinel_t && in_end, out_t & out_it, target_t & target)
{
    read_one(in_cur, in_end, out_it, seqan3::detail::equals_char<'o'>{});
    auto res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 1)), "h");
    EXPECT_EQ(*in_cur, 'e');

    read_one(in_cur, in_end, out_it, seqan3::detail::equals_char<'e'>{});
    res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 1)), "h");
    EXPECT_EQ(*in_cur, 'l');

    read_one(in_cur, in_end, out_it);
    res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 2)), "hl");
    EXPECT_EQ(*in_cur, 'l');
}

TYPED_TEST(tokenization_test, read_one)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;
    using target_type   = typename std::tuple_element<1, TypeParam>::type;
    using target_type_2 = std::conditional_t<std::is_same_v<target_type, std::ostringstream>,
                                             std::ostream, target_type>;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        target_type out{};
        auto it = test_tokenization_get_iterator(in);
        auto o_iter = make_preferred_output_iterator(out);

        impl_tokenization_test_read_one(it, test_tokenization_get_end(in), o_iter, out);
    }

    if constexpr (std::is_same_v<target_type_2, std::ostream>)
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        std::string storage;
        storage.resize(11);
        io_test_small_stream_buffer buf(storage.data(), storage.data() + storage.size());
        target_type_2 out{&buf};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_one(r_beg, r_end, o_iter, storage);
    }
    else
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        target_type_2 out{};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_one(r_beg, r_end, o_iter, out);
    }
}

template <typename in_t, typename sentinel_t, typename out_t, typename target_t>
void impl_tokenization_test_read_line(in_t & in_cur, sentinel_t && in_end, out_t & out_it, target_t & target)
{
    read_line(in_cur, in_end, out_it, seqan3::detail::equals_char<'o'>{});
    auto res = test_tokenization_get_data(target);
    EXPECT_EQ(static_cast<std::string>(res | ranges::view::slice(0, 9)), "hell_wrld");
    EXPECT_EQ(*in_cur, 't');

    read_line(in_cur, in_end, out_it);
    res = test_tokenization_get_data(target);
    EXPECT_EQ(res, "hell_wrldtest");
}

TYPED_TEST(tokenization_test, read_line)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;
    using target_type   = typename std::tuple_element<1, TypeParam>::type;
    using target_type_2 = std::conditional_t<std::is_same_v<target_type, std::ostringstream>,
                                             std::ostream, target_type>;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        target_type out{};
        auto it = test_tokenization_get_iterator(in);
        auto o_iter = make_preferred_output_iterator(out);

        impl_tokenization_test_read_line(it, test_tokenization_get_end(in), o_iter, out);
    }

    if constexpr (std::is_same_v<target_type_2, std::ostream>)
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        std::string storage;
        storage.resize(13);
        io_test_small_stream_buffer buf(storage.data(), storage.data() + storage.size());
        target_type_2 out{&buf};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_line(r_beg, r_end, o_iter, storage);
    }
    else
    {  // istream interface.
        source_type in{"hello_world\n\rtest"};
        target_type_2 out{};

        auto [r_beg, r_end] = make_preferred_input_iterator_range(in);
        auto o_iter = make_preferred_output_iterator(out);
        impl_tokenization_test_read_line(r_beg, r_end, o_iter, out);
    }
}

TYPED_TEST(tokenization_test, ignore_until)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        auto it = test_tokenization_get_iterator(in);
        read_until(it, test_tokenization_get_end(in), std::ignore, seqan3::detail::equals_char<'_'>{});
        EXPECT_EQ(*it, '_');
        read_until(it, test_tokenization_get_end(in), std::ignore, seqan3::detail::is_newline{});
        EXPECT_EQ(*it, '\n');
    }

    { // chunked
        source_type in{"hello_world\n\rtest"};
        auto [it, it_end] = make_preferred_input_iterator_range(in);
        read_until(it, it_end, std::ignore, seqan3::detail::equals_char<'_'>{});
        EXPECT_EQ(*it, '_');
        read_until(it, it_end, std::ignore, seqan3::detail::is_newline{});
        EXPECT_EQ(*it, '\n');
    }
}

TYPED_TEST(tokenization_test, ignore_n)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        auto it = test_tokenization_get_iterator(in);
        read_n(it, test_tokenization_get_end(in), std::ignore, 5);
        EXPECT_EQ(*it, '_');
        read_n(it, test_tokenization_get_end(in), std::ignore, 6);
        EXPECT_EQ(*it, '\n');
    }

    { // chunked
        source_type in{"hello_world\n\rtest"};
        auto [it, it_end] = make_preferred_input_iterator_range(in);
        read_n(it, it_end, std::ignore, 5);
        EXPECT_EQ(*it, '_');
        read_n(it, it_end, std::ignore, 6);
        EXPECT_EQ(*it, '\n');
    }
}

TYPED_TEST(tokenization_test, ignore_one)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        auto it = test_tokenization_get_iterator(in);
        read_one(it, test_tokenization_get_end(in), std::ignore);
        EXPECT_EQ(*it, 'e');
        read_one(it, test_tokenization_get_end(in), std::ignore);
        EXPECT_EQ(*it, 'l');
    }

    { // chunked
        source_type in{"hello_world\n\rtest"};
        auto [it, it_end] = make_preferred_input_iterator_range(in);
        read_one(it, it_end, std::ignore);
        EXPECT_EQ(*it, 'e');
        read_one(it, it_end, std::ignore);
        EXPECT_EQ(*it, 'l');
    }
}

TYPED_TEST(tokenization_test, ignore_line)
{
    using namespace seqan3::detail;

    using source_type   = typename std::tuple_element<0, TypeParam>::type;

    {  // non-chunked
        source_type in{"hello_world\n\rtest"};
        auto it = test_tokenization_get_iterator(in);
        read_line(it, test_tokenization_get_end(in), std::ignore);
        EXPECT_EQ(*it, 't');
    }

    { // chunked
        source_type in{"hello_world\n\rtest"};
        auto [it, it_end] = make_preferred_input_iterator_range(in);
        read_line(it, it_end, std::ignore);
        EXPECT_EQ(*it, 't');
    }
}
