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

#include <iostream>
#include <cstring>
#include <vector>

#ifndef SEQAN3_HAS_ZLIB
#error "This file cannot be used when building without ZLIB-support."
#endif

#include <zlib.h>

namespace seqan3::contrib
{

// Default gzip buffer size, change this to suite your needs.
const size_t GZ_INPUT_DEFAULT_BUFFER_SIZE = 921600;

// --------------------------------------------------------------------------
// Class basic_gz_istreambuf
// --------------------------------------------------------------------------
// A stream decorator that takes compressed input and unzips it to a istream.
// The class wraps up the deflate method of the zlib library 1.1.4 https://www.zlib.net

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_gz_istreambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_istream<Elem, Tr> &              istream_reference;
    typedef ElemA                                       char_allocator_type;
    typedef ByteT                                       byte_type;
    typedef ByteAT                                      byte_allocator_type;
    typedef byte_type *                                 byte_buffer_type;
    typedef Tr                                          traits_type;
    typedef typename Tr::char_type                      char_type;
    typedef typename Tr::int_type                       int_type;
    typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
    typedef std::vector<char_type, char_allocator_type> char_vector_type;

    // Construct a unzip stream
    // More info on the following parameters can be found in the zlib documentation.
    basic_gz_istreambuf(istream_reference istream_,
                          size_t window_size_,
                          size_t read_buffer_size_,
                          size_t input_buffer_size_);

    ~basic_gz_istreambuf();

    int_type underflow();

    // returns the compressed input istream
    istream_reference get_istream()  { return m_istream; }
    // returns the zlib stream structure
    z_stream & get_zip_stream()      { return m_zip_stream; }

private:
    void put_back_from_zip_stream();
    std::streamsize unzip_from_stream(char_type *, std::streamsize);
    size_t fill_input_buffer();

    istream_reference m_istream;
    z_stream m_zip_stream;
    int m_err;
    byte_vector_type m_input_buffer;
    char_vector_type m_buffer;
};

// --------------------------------------------------------------------------
// Class basic_gz_istreambuf implementation
// --------------------------------------------------------------------------

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::basic_gz_istreambuf(
    istream_reference istream_,
    size_t window_size_,
    size_t read_buffer_size_,
    size_t input_buffer_size_
    ) :
    m_istream(istream_),
    m_input_buffer(input_buffer_size_),
    m_buffer(read_buffer_size_)
{
    // setting zalloc, zfree and opaque
    m_zip_stream.zalloc = (alloc_func)0;
    m_zip_stream.zfree = (free_func)0;

    m_zip_stream.next_in = NULL;
    m_zip_stream.avail_in = 0;
    m_zip_stream.avail_out = 0;
    m_zip_stream.next_out = NULL;

    m_err = inflateInit2(&m_zip_stream, static_cast<int>(window_size_));

    this->setg(&(m_buffer[0]) + 4,  // beginning of putback area
               &(m_buffer[0]) + 4,  // read position
               &(m_buffer[0]) + 4); // end position
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::~basic_gz_istreambuf()
{
    inflateEnd(&m_zip_stream);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
typename basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type
basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::underflow()
{
    if (this->gptr() && (this->gptr() < this->egptr()))
        return *reinterpret_cast<unsigned char *>(this->gptr());

    int n_putback = static_cast<int>(this->gptr() - this->eback());
    if (n_putback > 4)
        n_putback = 4;

    std::memmove(&(m_buffer[0]) + (4 - n_putback), this->gptr() - n_putback, n_putback * sizeof(char_type));

    int num = unzip_from_stream(&(m_buffer[0]) + 4,
                                static_cast<std::streamsize>((m_buffer.size() - 4) * sizeof(char_type)));

    if (num <= 0)     // ERROR or EOF
        return traits_type::eof();

    // reset buffer pointers
    this->setg(&(m_buffer[0]) + (4 - n_putback),         // beginning of putback area
               &(m_buffer[0]) + 4,                       // read position
               &(m_buffer[0]) + 4 + num);                // end of buffer

    // return next character
    return *reinterpret_cast<unsigned char *>(this->gptr());
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::unzip_from_stream(
    char_type * buffer_,
    std::streamsize buffer_size_)
{
    m_zip_stream.next_out = (byte_buffer_type)buffer_;
    m_zip_stream.avail_out = static_cast<uInt>(buffer_size_ * sizeof(char_type));
    size_t count = m_zip_stream.avail_in;

    do
    {
        if (m_zip_stream.avail_in == 0)
            count = fill_input_buffer();

        if (m_zip_stream.avail_in)
            m_err = inflate(&m_zip_stream, Z_SYNC_FLUSH);

        if (m_err == Z_STREAM_END)
            inflateReset(&m_zip_stream);
        else if (m_err < 0)
            break;
    }
    while (m_zip_stream.avail_out > 0 && count > 0);

    std::streamsize n_read = buffer_size_ - m_zip_stream.avail_out / sizeof(char_type);

    // check if it is the end
    if (m_zip_stream.avail_out > 0 && m_err == Z_STREAM_END)
        put_back_from_zip_stream();

    return n_read;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
size_t basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::fill_input_buffer()
{
    m_zip_stream.next_in = &(m_input_buffer[0]);
    m_istream.read((char_type *)(&(m_input_buffer[0])),
                   static_cast<std::streamsize>(m_input_buffer.size() / sizeof(char_type)));
    return m_zip_stream.avail_in = m_istream.gcount() * sizeof(char_type);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
void basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT>::put_back_from_zip_stream()
{
    if (m_zip_stream.avail_in == 0)
        return;

    m_istream.clear(std::ios::goodbit);
    m_istream.seekg(-static_cast<int>(m_zip_stream.avail_in), std::ios_base::cur);

    m_zip_stream.avail_in = 0;
}

// --------------------------------------------------------------------------
// Class basic_gz_istreambase
// --------------------------------------------------------------------------
// Base class for unzip istreams
// Contains a basic_gz_istreambuf.

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_gz_istreambase :
    virtual public std::basic_ios<Elem, Tr>
{
public:
    typedef std::basic_istream<Elem, Tr> &                        istream_reference;
    typedef basic_gz_istreambuf<Elem, Tr, ElemA, ByteT, ByteAT> unzip_streambuf_type;

    basic_gz_istreambase(istream_reference ostream_,
                          size_t window_size_,
                          size_t read_buffer_size_,
                          size_t input_buffer_size_) :
        m_buf(ostream_, window_size_, read_buffer_size_, input_buffer_size_)
    {
        this->init(&m_buf);
    }

    // returns the underlying unzip istream object
    unzip_streambuf_type * rdbuf() { return &m_buf; }

private:
    unzip_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_gz_istream
// --------------------------------------------------------------------------
// A zipper istream
//
// This class is a istream decorator that behaves 'almost' like any other ostream.
// At construction, it takes any istream that shall be used to input of the compressed data.
//
// Simlpe example:
//
// // create a stream on zip string
// istringstream istringstream_( ostringstream_.str());
// // create unzipper istream
// zip_istream unzipper( istringstream_);
// // read and unzip
// unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_gz_istream :
    public basic_gz_istreambase<Elem, Tr, ElemA, ByteT, ByteAT>,
    public std::basic_istream<Elem, Tr>
{
public:
    typedef basic_gz_istreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_istreambase_type;
    typedef std::basic_istream<Elem, Tr>                          istream_type;
    typedef istream_type &                                        istream_reference;
    typedef ByteT                                                 byte_type;
    typedef Tr                                                    traits_type;

    // Construct a unzipper stream
    //
    // istream_ input buffer
    // window_size_
    // read_buffer_size_
    // input_buffer_size_

    basic_gz_istream(istream_reference istream_,
                      size_t window_size_ = 31, // 15 (size) + 16 (gzip header)
                      size_t read_buffer_size_ = GZ_INPUT_DEFAULT_BUFFER_SIZE,
                      size_t input_buffer_size_ = GZ_INPUT_DEFAULT_BUFFER_SIZE) :
        zip_istreambase_type(istream_, window_size_, read_buffer_size_, input_buffer_size_),
        istream_type(this->rdbuf())
    {}

#ifdef _WIN32
private:
    void _Add_vtordisp1() {}  // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() {}  // Required to avoid VC++ warning C4250
#endif
};

// ===========================================================================
// Typedefs
// ===========================================================================

// A typedef for basic_gz_istream<char>
typedef basic_gz_istream<char>     gz_istream;
// A typedef for basic_gz_istream<wchart>
typedef basic_gz_istream<wchar_t>  gz_wistream;

} // namespace seqan3::contrib
