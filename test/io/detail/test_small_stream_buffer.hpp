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

#include <array>
#include <streambuf>

class io_test_small_stream_buffer : public std::streambuf
{
public:

    using traits_type = typename std::streambuf::traits_type;

    io_test_small_stream_buffer(char * _begin, char * _end) : data_begin(_begin), data_end(_end)
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
