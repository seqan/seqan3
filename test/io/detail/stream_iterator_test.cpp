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

// TODO(rrahn): Make typed test for different containers?
#include <gtest/gtest.h>

#include <cstring>
#include <cassert>

#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <array>
#include <streambuf>

#include <seqan3/io/detail/stream_iterator.hpp>

namespace seqan3::detail::test_stream_iterator
{

class small_stream_buffer : public std::streambuf
{
public:

    using traits_type = typename std::streambuf::traits_type;

    small_stream_buffer(char * _begin, char * _end) : data_begin(_begin), data_end(_end)
    {
        assert(_end - _begin >= static_cast<std::ptrdiff_t>(BUFFER_SIZE));
        setg(_begin, _begin, _begin + BUFFER_SIZE);
        setp(_begin, _begin + BUFFER_SIZE);
    }

protected:

    typename traits_type::int_type underflow()
    {
        using std::size;

        // Valid position, return current char.
        if (gptr() < egptr())
            return traits_type::to_int_type(*gptr());

        // Reached end of stream.
        if (gptr() == data_end)
            return traits_type::eof();

        // Store data last chars in putback buffer.
        std::memmove(&put_back_buffer[0], egptr() - size(put_back_buffer), size(put_back_buffer));

        // Update buffer.
        setg(gptr(), gptr(), gptr() + std::min(static_cast<std::ptrdiff_t>(BUFFER_SIZE), data_end - gptr()));
        return traits_type::to_int_type(*gptr());
    }

    typename traits_type::int_type overflow(typename traits_type::int_type ch = traits_type::eof())
    {
        // reset the put area.
        if (pptr() == epptr())
            setp(pptr(), pptr() + std::min(static_cast<std::ptrdiff_t>(BUFFER_SIZE), data_end - pptr()));

        if (pptr() == data_end || traits_type::eq_int_type(ch, traits_type::eof()))
            return traits_type::eof();

        *pptr() = ch;
        return ch;
    }
private:

    static const size_t BUFFER_SIZE = 3;

    char * data_begin{nullptr};
    char * data_end{nullptr};

    std::array<char, 1>  put_back_buffer;
};

}  // seqan3::detail::test_stream_iterator

using namespace seqan3;

// What we actually want is a direction_iterator over stl containers and
// And we want to make them applicable to streams
TEST(stream_iterator, chunk_istream_iterator_construction)
{
    using namespace seqan3::detail;

    std::stringstream in{"acg"};

    EXPECT_NO_THROW(chunk_istream_iterator<std::stringstream>{in});
    EXPECT_NO_THROW(chunk_istream_iterator<std::stringstream>{in.rdbuf()});
    EXPECT_NO_THROW(chunk_istream_iterator<std::stringstream>{});
}

TEST(stream_iterator, chunk_istream_iterator_dereference)
{
    using namespace seqan3::detail;

    std::stringstream in{"acgtgatagctacgacgatcg"};

    chunk_istream_iterator<std::stringstream> it{in};

    EXPECT_EQ(*it, 'a');
}

TEST(stream_iterator, chunk_istream_iterator_increment)
{
    using namespace seqan3::detail;

    std::stringstream in{"acgtgatagctacgacgatcg"};

    chunk_istream_iterator<std::stringstream> it{in};

    ++it;
    EXPECT_EQ(*it, 'c');
    it++;
    EXPECT_EQ(*it, 'g');
}

TEST(stream_iterator, chunk_istream_iterator_get_chunk)
{
    using namespace seqan3::detail;

    {
        std::string in{"acgtgatagctacgacgatcg"};
        test_stream_iterator::small_stream_buffer buf(in.data(), in.data() + in.size());

        chunk_istream_iterator<std::iostream> it{&buf};

        auto rng = it.get_chunk();

        EXPECT_EQ(*rng.begin(), 'a');
        EXPECT_NE(rng.begin(), rng.end());
        EXPECT_EQ(rng.end() - rng.begin(), 3);
    }

    {  // empty stream
        std::stringstream in{""};

        chunk_istream_iterator<std::stringstream> it{in};

        auto rng = it.get_chunk();

        EXPECT_EQ(rng.begin(), rng.end());
        EXPECT_EQ(rng.end() - rng.begin(), 0);
    }
}

TEST(stream_iterator, chunk_istream_iterator_next_chunk)
{
    using namespace seqan3::detail;

    {  // empty stream
        std::stringstream in{""};

        chunk_istream_iterator<std::stringstream> it{in};

        it.next_chunk();
        auto rng = it.get_chunk();

        EXPECT_EQ(rng.begin(), rng.end());
        EXPECT_EQ(rng.end() - rng.begin(), 0);
    }

    {
        std::string in{"acgtgata"};
        test_stream_iterator::small_stream_buffer buf(in.data(), in.data() + in.size());

        chunk_istream_iterator<std::iostream> it{&buf};

        // Is not doing anything, because we are not at the end of the stream.
        it.next_chunk();
        auto rng = it.get_chunk();

        EXPECT_EQ(*rng.begin(), in[0]);
        EXPECT_EQ(rng.end() - rng.begin(), 3);

        advance(it, 3);

        it.next_chunk();
        rng = it.get_chunk();

        EXPECT_EQ(*rng.begin(), in[3]);
        EXPECT_EQ(rng.end() - rng.begin(), 3);

        advance(it, 3);

        it.next_chunk();
        rng = it.get_chunk();
        EXPECT_EQ(*rng.begin(), in[6]);
        EXPECT_EQ(rng.end() - rng.begin(), 2);
    }
}

TEST(stream_iterator, chunk_istream_iterator_advance_chunk)
{
    using namespace seqan3::detail;

    using std::advance;

    std::string in{"acgtgata"};
    test_stream_iterator::small_stream_buffer buf(in.data(), in.data() + in.size());

    chunk_istream_iterator<std::iostream> it{&buf};

    it.advance_chunk(2);
    auto rng = it.get_chunk();

    EXPECT_EQ(*rng.begin(), in[2]);
    EXPECT_EQ(rng.end() - rng.begin(), 1);

    it.advance_chunk(1);
    it.next_chunk();
    it.advance_chunk(1);
    rng = it.get_chunk();

    EXPECT_EQ(*rng.begin(), in[4]);
    EXPECT_EQ(rng.end() - rng.begin(), 2);
}

TEST(stream_iterator, chunk_ostream_iterator_construction)
{
    using namespace seqan3::detail;

    std::ostringstream out{};

    EXPECT_NO_THROW(chunk_ostream_iterator<std::ostringstream>{out});
    EXPECT_NO_THROW(chunk_ostream_iterator<std::ostringstream>{out.rdbuf()});
    EXPECT_NO_THROW(chunk_ostream_iterator<std::ostringstream>{});
}

TEST(stream_iterator, chunk_ostream_iterator_increment)
{
    using namespace seqan3::detail;

    std::ostringstream out{};
    chunk_ostream_iterator<std::ostringstream> it{out};
    EXPECT_NO_THROW(it++);
    EXPECT_NO_THROW(++it);
    EXPECT_EQ(out.str(), "");
}

TEST(stream_iterator, chunk_ostream_iterator_dereference)
{
    using namespace seqan3::detail;

    std::ostringstream out{};
    chunk_ostream_iterator<std::ostringstream> it{out};
    EXPECT_NO_THROW(*it);
    EXPECT_TRUE((std::is_same_v<decltype(*it), chunk_ostream_iterator<std::ostringstream>&>));
}

TEST(stream_iterator, chunk_ostream_iterator_assign)
{
    using namespace seqan3::detail;

    std::ostringstream out{};

    chunk_ostream_iterator<std::ostringstream> it{out};

    auto& ret = it = 'v';

    EXPECT_TRUE((std::is_same_v<decltype(it = 'a'), chunk_ostream_iterator<std::ostringstream>&>));
    ret = 'a';
    EXPECT_EQ(out.str(), "va");
}

TEST(stream_iterator, chunk_ostream_iterator_get_chunk)
{
    using namespace seqan3::detail;

    std::string out;
    out.resize(10);
    test_stream_iterator::small_stream_buffer buf(out.data(), out.data() + out.size());

    chunk_ostream_iterator<std::iostream> it{&buf};

    auto rng = it.get_chunk();

    *rng.begin() = 'a';
    *(rng.begin() + 1) = 'c';
    *(rng.begin() + 2) = 'g';

    EXPECT_EQ(out[0], 'a');
    EXPECT_EQ(out[1], 'c');
    EXPECT_EQ(out[2], 'g');
    EXPECT_EQ(rng.end() - rng.begin(), 3);
}

TEST(stream_iterator, chunk_ostream_iterator_next_chunk)
{
    using namespace seqan3::detail;

    std::string out;
    out.resize(10);
    test_stream_iterator::small_stream_buffer buf(out.data(), out.data() + out.size());

    chunk_ostream_iterator<std::iostream> it{&buf};

    it = 'a';
    it = 'b';

    EXPECT_EQ(out[0], 'a');
    EXPECT_EQ(out[1], 'b');

    it.next_chunk();

    auto rng = it.get_chunk();

    EXPECT_EQ(rng.begin(), out.data() + 2);
    EXPECT_EQ(rng.end(), out.data() + 3);

    it = 'c';
    EXPECT_EQ(out[2], 'c');

    it.next_chunk();

    rng = it.get_chunk();

    EXPECT_EQ(rng.begin(), out.data() + 3);
    EXPECT_EQ(rng.end(), out.data() + 6);
}

TEST(stream_iterator, chunk_ostream_iterator_advance_chunk)
{
    using namespace seqan3::detail;

    std::string out;
    out.resize(10);
    test_stream_iterator::small_stream_buffer buf(out.data(), out.data() + out.size());

    chunk_ostream_iterator<std::iostream> it{&buf};

    it.advance_chunk(2);
    it.next_chunk();
    it = 'a';
    EXPECT_EQ(out[2], 'a');

    it.next_chunk();
    it = 'a';
    it = 'c';
    EXPECT_EQ(out[3], 'a');
    EXPECT_EQ(out[4], 'c');

    it.advance_chunk(1);
    it.next_chunk();
    it = 'c';
    EXPECT_EQ(out[6], 'c');
}

TEST(stream_iterator, input_iterator_stream)
{
    using namespace seqan3::detail;

    std::stringstream in{"acg"};
    auto [r_beg, r_end] = input_iterator(in);
    EXPECT_EQ(*r_beg, 'a');
    EXPECT_EQ(r_end, chunk_istream_iterator<std::stringstream>{});
}

TEST(stream_iterator, output_iterator_stream)
{
    using namespace seqan3::detail;

    std::stringstream out{};
    out << "acg";

    *output_iterator(out) = 't';
    EXPECT_EQ(out.str()[3], 't');
    EXPECT_EQ(out.str().size(), 4u);
}
