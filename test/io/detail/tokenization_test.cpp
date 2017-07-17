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

#include <range/v3/utility/iterator.hpp>

#include <seqan3/io/detail/tokenization.hpp>

using namespace seqan3;

// What we actually want is a direction_iterator over stl containers and
// And we want to make them applicable to streams

// For example our own stream iterator can do this, if we pass our own stream_adapter.
// We want them chunkable if the container allows for it.

TEST(tokenization_test, write_container_to_container)
{
    using namespace seqan3::detail;

    using std::size;
    using std::begin;

    std::string out{};
    std::string in{"hello world"};

    {  // standrd iterator + output_iterator
        write(begin(in), size(in), output_iterator(out));
        EXPECT_EQ(out, in);
    }

    {  // standrd iterator + output
        out.clear();
        write(begin(in), size(in), out);
        EXPECT_EQ(out, in);
    }

    {  // input_iterator + output_iterator
        out.clear();
        write(std::get<0>(input_iterator(in)), size(in), output_iterator(out));
        EXPECT_EQ(out, in);
    }

    {  // container interface
        out.clear();
        write(in, out);
        EXPECT_EQ(out, in);
    }
}

TEST(tokenization_test, write_container_to_stream)
{
    using namespace seqan3::detail;

    std::string in{"hello world"};

    {  // Container input interface.
        std::ostringstream out{};
        write(in, out);
        EXPECT_EQ(out.str(), in);
    }

    {  // standard iterator interface.
        std::ostringstream out{};
        write(begin(in), size(in), out);
        EXPECT_EQ(out.str(), in);
    }

    {  // standard iterator interface + output_iterator.
        std::ostringstream out{};
        write(begin(in), size(in), output_iterator(out));
        EXPECT_EQ(out.str(), in);
    }

    {  // input_iterator interface + output_iterator.
        std::ostringstream out{};
        write(std::get<0>(input_iterator(in)), size(in), output_iterator(out));
        EXPECT_EQ(out.str(), in);
    }
}

TEST(tokenization_test, write_stream_to_container)
{
    using namespace seqan3::detail;

    {  // istream interface.
        std::istringstream in{"hello_world"};
        std::string out;
        std::istream_iterator<char> it{in};
        write(it, 11, out);
        EXPECT_EQ(out, in.str());
    }

    {  // input_iterator interface.
        std::istringstream in{"hello_world"};
        std::string out;
        write(std::get<0>(input_iterator(in)), 11, out);
        EXPECT_EQ(out, in.str());
    }
}

TEST(tokenization_test, write_stream_to_stream)
{
    using namespace seqan3::detail;

    {  // istream interface.
        std::istringstream in{"hello_world"};
        std::ostringstream out;
        std::istream_iterator<char> i_it{in};

        write(i_it, 11, out);
        EXPECT_EQ(out.str(), in.str());
    }

    {  // istream interface.
        std::istringstream in{"hello_world"};
        std::ostringstream out;
        std::istream_iterator<char> i_it{in};

        // ostream iterator does not implement the value_type metafunction -> evaluates to void?
        ranges::v3::ostream_iterator<char> o_it{out};

        write(i_it, 11, o_it);
        EXPECT_EQ(out.str(), in.str());
    }

    {  // input_iterator interface.
        std::istringstream in{"hello_world"};
        std::ostringstream out;
        write(std::get<0>(input_iterator(in)), 11, output_iterator(out));
        EXPECT_EQ(out.str(), in.str());
    }
}

TEST(tokenization_test, write_container_to_array)
{

}

TEST(tokenization_test, write_array_to_stream)
{

}

TEST(tokenization_test, write_stream_to_array)
{

}

TEST(tokenization_test, write_array_to_array)
{

}
