// bzip2stream Library License:
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
// Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
// Altered bzip2_stream header
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>

#pragma once

#ifndef SEQAN3_HAS_BZIP2
#error "This file cannot be used when building without BZIP2-support."
#endif

#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

#define BZ_NO_STDIO
#include <bzlib.h>

namespace seqan3::contrib
{

// --------------------------------------------------------------------------
// Class basic_bz2_istreambuf
// --------------------------------------------------------------------------

const size_t BZ2_INPUT_DEFAULT_BUFFER_SIZE = 4096;

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bz2_istreambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef ElemA char_allocator_type;
    typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
    typedef byte_type* byte_buffer_type;
    typedef typename Tr::char_type char_type;
    typedef typename Tr::int_type int_type;
    typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
    typedef std::vector<char_type, char_allocator_type > char_vector_type;

    basic_bz2_istreambuf(
        istream_reference istream_,
        size_t verbosity_,
        bool small_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
        );

    ~basic_bz2_istreambuf();

    int_type underflow();

    istream_reference get_istream()        {    return m_istream;};
    bz_stream& get_bzip2_stream()        {    return m_bzip2_stream;};
    int get_zerr() const                {    return m_err;};
private:
    std::streamsize unbzip2_from_stream( char_type*, std::streamsize);
    void put_back_from_bzip2_stream();
    size_t fill_input_buffer();

    istream_reference m_istream;
    bz_stream m_bzip2_stream;
    int m_err;
    byte_vector_type m_input_buffer;
    char_vector_type m_buffer;
};

// --------------------------------------------------------------------------
// Class basic_bz2_istreambuf implementation
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
basic_bz2_istreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::basic_bz2_istreambuf(
        istream_reference istream_,
        size_t verbosity_,
        bool small_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
)
:
    m_istream(istream_),
    m_input_buffer(input_buffer_size_),
    m_buffer(read_buffer_size_)
{
    // setting zalloc, zfree and opaque
    m_bzip2_stream.bzalloc=NULL;
    m_bzip2_stream.bzfree=NULL;

    m_bzip2_stream.next_in=NULL;
    m_bzip2_stream.avail_in=0;
    m_bzip2_stream.avail_out=0;
    m_bzip2_stream.next_out=NULL;


    m_err=BZ2_bzDecompressInit (
        &m_bzip2_stream,
        std::min(4, static_cast<int>(verbosity_)),
        static_cast<int>(small_)
    );

    this->setg(
        &(m_buffer[0])+4,     // beginning of putback area
        &(m_buffer[0])+4,     // read position
        &(m_buffer[0])+4);    // end position
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
size_t basic_bz2_istreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::fill_input_buffer()
{
    m_bzip2_stream.next_in=&(m_input_buffer[0]);
    m_istream.read(
        (char_type*)(&(m_input_buffer[0])),
        static_cast<std::streamsize>(m_input_buffer.size()/sizeof(char_type))
        );
    return m_bzip2_stream.avail_in=m_istream.gcount()*sizeof(char_type);
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
void basic_bz2_istreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::put_back_from_bzip2_stream()
{
    if (m_bzip2_stream.avail_in==0)
        return;

    m_istream.clear( std::ios::goodbit );
    m_istream.seekg(
        -static_cast<int>(m_bzip2_stream.avail_in),
        std::ios_base::cur
        );

    m_bzip2_stream.avail_in=0;
}


template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
basic_bz2_istreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::~basic_bz2_istreambuf()
{
    BZ2_bzDecompressEnd(&m_bzip2_stream);
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
typename basic_bz2_istreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::int_type
        basic_bz2_istreambuf<
            Elem,Tr,ElemA,ByteT,ByteAT
            >::underflow()
{
    if ( this->gptr() && ( this->gptr() < this->egptr()))
        return * reinterpret_cast<unsigned char *>( this->gptr());

    int n_putback = static_cast<int>(this->gptr() - this->eback());
    if ( n_putback > 4)
        n_putback = 4;
    std::memmove(
        &(m_buffer[0]) + (4 - n_putback),
        this->gptr() - n_putback,
        n_putback*sizeof(char_type)
        );

    int num = unbzip2_from_stream(
        &(m_buffer[0])+4,
        static_cast<std::streamsize>((m_buffer.size()-4)*sizeof(char_type))
        );
    if (num <= 0) // ERROR or EOF
        return EOF;

    // reset buffer pointers
    this->setg(
            &(m_buffer[0]) + (4 - n_putback),   // beginning of putback area
            &(m_buffer[0]) + 4,                 // read position
            &(m_buffer[0]) + 4 + num);          // end of buffer

        // return next character
        return* reinterpret_cast<unsigned char *>( this->gptr());
    }


template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
std::streamsize basic_bz2_istreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::unbzip2_from_stream(
        char_type* buffer_,
        std::streamsize buffer_size_
        )
{
    m_bzip2_stream.next_out=(byte_buffer_type)buffer_;
    m_bzip2_stream.avail_out=buffer_size_*sizeof(char_type);
    size_t count =m_bzip2_stream.avail_in;

    do
    {
        if (m_bzip2_stream.avail_in==0)
            count=fill_input_buffer();

        if (m_bzip2_stream.avail_in)
        {
            m_err = BZ2_bzDecompress( &m_bzip2_stream );
        }
    } while (m_err==BZ_OK && m_bzip2_stream.avail_out != 0 && count != 0);

    if (m_err == BZ_STREAM_END)
        put_back_from_bzip2_stream();

    return buffer_size_ - m_bzip2_stream.avail_out/sizeof(char_type);
}

// --------------------------------------------------------------------------
// Class basic_bz2_istreambase
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bz2_istreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_istream<Elem, Tr>& istream_reference;
    typedef basic_bz2_istreambuf<
        Elem,Tr,ElemA,ByteT,ByteAT> unbzip2_streambuf_type;

    basic_bz2_istreambase(
        istream_reference ostream_,
        size_t verbosity_,
        bool small_,
        size_t read_buffer_size_,
        size_t input_buffer_size_
        )
        : m_buf(
            ostream_,
            verbosity_,
            small_,
            read_buffer_size_,
            input_buffer_size_
            )
    {
        this->init(&m_buf );
    };

    unbzip2_streambuf_type* rdbuf() { return &m_buf; };

private:
    unbzip2_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_bz2_istream
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bz2_istream :
    public basic_bz2_istreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_istream<Elem,Tr>
{
public:
    typedef basic_bz2_istreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bzip2_istreambase_type;
    typedef std::basic_istream<Elem,Tr> istream_type;
    typedef istream_type& istream_reference;
    typedef unsigned char byte_type;

    basic_bz2_istream(
        istream_reference istream_,
        size_t verbosity_ = 0,
        bool small_ = false,
        size_t read_buffer_size_ = BZ2_INPUT_DEFAULT_BUFFER_SIZE,
        size_t input_buffer_size_ = BZ2_INPUT_DEFAULT_BUFFER_SIZE
        )
      :
        bzip2_istreambase_type(istream_,verbosity_, small_, read_buffer_size_, input_buffer_size_),
        istream_type(bzip2_istreambase_type::rdbuf())
    {};
#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

// --------------------------------------------------------------------------
// typedefs
// --------------------------------------------------------------------------

typedef basic_bz2_istream<char> bz2_istream;
typedef basic_bz2_istream<wchar_t> bz2_wistream;

} // namespace seqan3::contrib
