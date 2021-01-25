// zipstream Library License:
// --------------------------
//
// The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.
//
// This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
//
// 3. This notice may not be removed or altered from any source distribution
//
// Altered zipstream library header
// Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>

#pragma once

#ifndef SEQAN3_HAS_ZLIB
#error "This file cannot be used when building without ZLIB-support."
#endif

#include <iostream>
#include <cstring>
#include <vector>

#include <zlib.h>

namespace seqan3::contrib
{

// Default gzip buffer size, change this to suite your needs.
const size_t GZ_OUTPUT_DEFAULT_BUFFER_SIZE = 921600;

// --------------------------------------------------------------------------
// Enum EStrategy
// --------------------------------------------------------------------------
// Compression strategy, see zlib doc.

enum EStrategy
{
    StrategyFiltered = 1,
    StrategyHuffmanOnly = 2,
    DefaultStrategy = 0
};

// --------------------------------------------------------------------------
// Class basic_gz_ostreambuf
// --------------------------------------------------------------------------
// A stream decorator that takes raw input and zips it to a ostream.
// The class wraps up the inflate method of the zlib library 1.1.4 https://www.zlib.net

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_gz_ostreambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr> &              ostream_reference;
    typedef ElemA                                       char_allocator_type;
    typedef ByteT                                       byte_type;
    typedef ByteAT                                      byte_allocator_type;
    typedef byte_type *                                 byte_buffer_type;
    typedef Tr                                          traits_type;
    typedef typename Tr::char_type                      char_type;
    typedef typename Tr::int_type                       int_type;
    typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
    typedef std::vector<char_type, char_allocator_type> char_vector_type;

    // Construct a zip stream
    // More info on the following parameters can be found in the zlib documentation.
    basic_gz_ostreambuf(ostream_reference ostream_,
                        size_t level_,
                        EStrategy strategy_,
                        size_t window_size_,
                        size_t memory_level_,
                        size_t buffer_size_);

    ~basic_gz_ostreambuf();

    int sync();
    int_type overflow(int_type c);

    // flushes the zip buffer and output buffer.
    // This method should be called at the end of the compression.
    // Calling flush multiple times, will lower the compression ratio.
    std::streamsize flush();

    // flushes the zip buffer and output buffer and finalize the zip stream
    // This method should be called at the end of the compression.
    std::streamsize flush_finalize();


private:
    bool zip_to_stream(char_type *, std::streamsize);
    size_t fill_input_buffer();
    // flush the zip buffer using a particular mode and flush output buffer
    std::streamsize flush(int flush_mode);

    ostream_reference m_ostream;
    z_stream m_zip_stream;
    int m_err;
    byte_vector_type m_output_buffer;
    char_vector_type m_buffer;
};

// --------------------------------------------------------------------------
// Class basic_gz_ostreambuf implementation
// --------------------------------------------------------------------------

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::basic_gz_ostreambuf(
    ostream_reference ostream_,
    size_t level_,
    EStrategy strategy_,
    size_t window_size_,
    size_t memory_level_,
    size_t buffer_size_
    ) :
    m_ostream(ostream_),
    m_output_buffer(buffer_size_, 0),
    m_buffer(buffer_size_, 0)
{
    m_zip_stream.zalloc = (alloc_func)0;
    m_zip_stream.zfree = (free_func)0;

    m_zip_stream.next_in = NULL;
    m_zip_stream.avail_in = 0;
    m_zip_stream.avail_out = 0;
    m_zip_stream.next_out = NULL;

    m_err = deflateInit2(
        &m_zip_stream,
        std::min(9, static_cast<int>(level_)),
        Z_DEFLATED,
        static_cast<int>(window_size_),
        std::min(9, static_cast<int>(memory_level_)),
        static_cast<int>(strategy_)
        );

    this->setp(&(m_buffer[0]), &(m_buffer[m_buffer.size() - 1]));
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::~basic_gz_ostreambuf()
{
    flush_finalize();
    m_ostream.flush();
    m_err = deflateEnd(&m_zip_stream);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
int basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::sync()
{
    if (this->pptr() && this->pptr() > this->pbase())
    {
        if (traits_type::eq_int_type(overflow(traits_type::eof()), traits_type::eof()))
            return -1;
    }

    return 0;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
typename basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type
basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::overflow(
    typename basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type c)
{
    int w = static_cast<int>(this->pptr() - this->pbase());

    if (!traits_type::eq_int_type(c, traits_type::eof()))
    {
        *this->pptr() = c;
        ++w;
    }

    if (zip_to_stream(this->pbase(), w))
    {
        this->setp(this->pbase(), this->epptr() - 1);
        return c;
    }
    else
    {
        return traits_type::eof();
    }
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
bool basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::zip_to_stream(
    typename basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::char_type * buffer_,
    std::streamsize buffer_size_)
{
    std::streamsize written_byte_size = 0, total_written_byte_size = 0;

    m_zip_stream.next_in = (byte_buffer_type)buffer_;
    m_zip_stream.avail_in = static_cast<uInt>(buffer_size_ * sizeof(char_type));
    m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size());
    m_zip_stream.next_out = &(m_output_buffer[0]);
    size_t remainder = 0;

    do
    {
        m_err = deflate(&m_zip_stream, 0);

        if (m_err == Z_OK  || m_err == Z_STREAM_END)
        {
            written_byte_size = static_cast<std::streamsize>(m_output_buffer.size()) - m_zip_stream.avail_out;
            total_written_byte_size += written_byte_size;

            // output buffer is full, dumping to ostream
            m_ostream.write((const char_type *) &(m_output_buffer[0]),
                            static_cast<std::streamsize>(written_byte_size / sizeof(char_type)));

            // checking if some bytes were not written.
            if ((remainder = written_byte_size % sizeof(char_type)) != 0)
            {
                // copy to the beginning of the stream
                std::memmove(&(m_output_buffer[0]),
                             &(m_output_buffer[written_byte_size - remainder]),
                             remainder);
            }

            m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size() - remainder);
            m_zip_stream.next_out = &m_output_buffer[remainder];
        }
    }
    while (m_zip_stream.avail_in != 0 && m_err == Z_OK);

    return m_err == Z_OK;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::flush(int flush_mode)
{
    int const buffer_size = static_cast<int>(this->pptr() - this->pbase()); // amount of data currently in buffer

    std::streamsize written_byte_size = 0, total_written_byte_size = 0;

    m_zip_stream.next_in = (byte_buffer_type) this->pbase();
    m_zip_stream.avail_in = static_cast<uInt>(buffer_size * sizeof(char_type));
    m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size());
    m_zip_stream.next_out = &(m_output_buffer[0]);
    size_t remainder = 0;

    do
    {
        m_err = deflate(&m_zip_stream, flush_mode);
        if (m_err == Z_OK || m_err == Z_STREAM_END)
        {
            written_byte_size = static_cast<std::streamsize>(m_output_buffer.size()) - m_zip_stream.avail_out;
            total_written_byte_size += written_byte_size;

            // output buffer is full, dumping to ostream
            m_ostream.write((const char_type *) &(m_output_buffer[0]),
                            static_cast<std::streamsize>(written_byte_size / sizeof(char_type) * sizeof(byte_type)));

            // checking if some bytes were not written.
            if ((remainder = written_byte_size % sizeof(char_type)) != 0)
            {
                // copy to the beginning of the stream
                std::memmove(&(m_output_buffer[0]),
                             &(m_output_buffer[written_byte_size - remainder]),
                             remainder);
            }

            m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size() - remainder);
            m_zip_stream.next_out = &m_output_buffer[remainder];
        }
    }
    while (m_err == Z_OK);

    m_ostream.flush();

    return total_written_byte_size;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::flush()
{
    return flush(Z_SYNC_FLUSH);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::flush_finalize()
{
    return flush(Z_FINISH);
}

// --------------------------------------------------------------------------
// Class basic_gz_ostreambase
// --------------------------------------------------------------------------
// Base class for zip ostreams.
// Contains a basic_gz_ostreambuf.

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_gz_ostreambase :
    virtual public std::basic_ios<Elem, Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr> &                      ostream_reference;
    typedef basic_gz_ostreambuf<Elem, Tr, ElemA, ByteT, ByteAT> zip_streambuf_type;

    // Construct a zip stream
    // More info on the following parameters can be found in the zlib documentation.
    basic_gz_ostreambase(ostream_reference ostream_,
                          size_t level_,
                          EStrategy strategy_,
                          size_t window_size_,
                          size_t memory_level_,
                          size_t buffer_size_) :
        m_buf(ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_)
    {
        this->init(&m_buf);
    }

    // returns the underlying zip ostream object
    zip_streambuf_type * rdbuf() { return &m_buf; }

private:
    zip_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_gz_ostream
// --------------------------------------------------------------------------
// A zipper ostream
//
// This class is a ostream decorator that behaves 'almost' like any other ostream.
// At construction, it takes any ostream that shall be used to output of the compressed data.
// When finished, you need to call the special method zflush or call the destructor
// to flush all the intermidiate streams.
//
// Example:
//
// // creating the target zip string, could be a fstream
// ostringstream ostringstream_;
// // creating the zip layer
// zip_ostream zipper(ostringstream_);
// // writing data
// zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// // zip ostream needs special flushing...
// zipper.zflush();

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_gz_ostream :
    public basic_gz_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT>,
    public std::basic_ostream<Elem, Tr>
{
public:
    typedef basic_gz_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_ostreambase_type;
    typedef std::basic_ostream<Elem, Tr>                          ostream_type;
    typedef ostream_type &                                        ostream_reference;

    // Constructs a zipper ostream decorator
    //
    // ostream_ ostream where the compressed output is written
    // is_gzip_ true if gzip header and footer have to be added
    // level_ level of compression 0, bad and fast, 9, good and slower,
    // strategy_ compression strategy
    // window_size_ see zlib doc
    // memory_level_ see zlib doc
    // buffer_size_ the buffer size used to zip data

    basic_gz_ostream(ostream_reference ostream_,
                      size_t level_ = Z_DEFAULT_COMPRESSION,
                      EStrategy strategy_ = DefaultStrategy,
                      size_t window_size_ = 31, // 15 (size) + 16 (gzip header)
                      size_t memory_level_ = 8,
                      size_t buffer_size_ = GZ_OUTPUT_DEFAULT_BUFFER_SIZE) :
        zip_ostreambase_type(ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_),
        ostream_type(this->rdbuf())
    {}

    ~basic_gz_ostream()
    {
        ostream_type::flush(); this->rdbuf()->flush_finalize();
    }

    // flush inner buffer and zipper buffer
    basic_gz_ostream<Elem, Tr> & flush()
    {
        ostream_type::flush(); this->rdbuf()->flush(); return *this;
    }

#ifdef _WIN32
private:
    void _Add_vtordisp1() {}  // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() {}  // Required to avoid VC++ warning C4250
#endif
};

// ===========================================================================
// Typedefs
// ===========================================================================

// A typedef for basic_gz_ostream<char>
typedef basic_gz_ostream<char>     gz_ostream;
// A typedef for basic_gz_ostream<wchar_t>
typedef basic_gz_ostream<wchar_t>  gz_wostream;

} // namespace seqan3::contrib
