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

#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>

#ifndef SEQAN3_HAS_BZIP2
#error "This file cannot be used when building without BZIP2-support."
#endif

#define BZ_NO_STDIO
#include <bzlib.h>

namespace seqan3::contrib
{

// --------------------------------------------------------------------------
// Class basic_bz2_ostreambuf
// --------------------------------------------------------------------------

const size_t BZ2_OUTPUT_DEFAULT_BUFFER_SIZE = 4096;

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bz2_ostreambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_streambuf< Elem, Tr > basic_streambuf_type;
    typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef ElemA char_allocator_type;
    typedef ByteT byte_type;
    typedef ByteAT byte_allocator_type;
    typedef byte_type* byte_buffer_type;
    typedef typename Tr::char_type char_type;
    typedef typename Tr::int_type int_type;
    typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
    typedef std::vector<char_type, char_allocator_type > char_vector_type;

    using basic_streambuf_type::epptr;
    using basic_streambuf_type::pbase;
    using basic_streambuf_type::pptr;

    basic_bz2_ostreambuf(
        ostream_reference ostream_,
        size_t block_size_100k_ ,
        size_t verbosity_ ,
        size_t work_factor_,
        size_t buffer_size_
        );

    ~basic_bz2_ostreambuf();

    int sync ();
    int_type overflow (int_type c);

    std::streamsize flush(int flush_mode);
    int get_zerr() const
    {    return m_err;};
    uint64_t get_in_size() const
    {
        return ((uint64_t)m_bzip2_stream.total_in_hi32 << 32)
                + m_bzip2_stream.total_in_lo32;
    }
    uint64_t get_out_size() const
    {
        return ((uint64_t)m_bzip2_stream.total_out_hi32 << 32)
                + m_bzip2_stream.total_out_lo32;
    }
private:
    bool bzip2_to_stream( char_type*, std::streamsize);
    size_t fill_input_buffer();

    ostream_reference m_ostream;
    bz_stream m_bzip2_stream;
    int m_err;
    byte_vector_type m_output_buffer;
    char_vector_type m_buffer;
};

// --------------------------------------------------------------------------
// Class basic_bz2_ostreambuf implementation
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
basic_bz2_ostreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >:: basic_bz2_ostreambuf(
    ostream_reference ostream_,
    size_t block_size_100k_,
    size_t verbosity_,
    size_t work_factor_,
    size_t buffer_size_
    )
:
    m_ostream(ostream_),
    m_output_buffer(buffer_size_,0),
    m_buffer(buffer_size_,0)
{
    m_bzip2_stream.bzalloc=NULL;
    m_bzip2_stream.bzfree=NULL;

    m_bzip2_stream.next_in=NULL;
    m_bzip2_stream.avail_in=0;
    m_bzip2_stream.avail_out=0;
    m_bzip2_stream.next_out=NULL;

    m_err=BZ2_bzCompressInit(
        &m_bzip2_stream,
        std::min( 9, static_cast<int>(block_size_100k_) ),
        std::min( 4, static_cast<int>(verbosity_) ),
        std::min( 250, static_cast<int>(work_factor_) )
        );

    this->setp( &(m_buffer[0]), &(m_buffer[m_buffer.size()-1]));
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
basic_bz2_ostreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::~basic_bz2_ostreambuf()
{
    flush(BZ_FINISH);
    m_ostream.flush();
    m_err=BZ2_bzCompressEnd(&m_bzip2_stream);
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
int basic_bz2_ostreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::sync ()
{
    if ( this->pptr() && this->pptr() > this->pbase())
    {
        int c = overflow( EOF);

        if ( c == EOF)
            return -1;
    }

    return 0;
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
typename basic_bz2_ostreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::int_type
        basic_bz2_ostreambuf<
            Elem,Tr,ElemA,ByteT,ByteAT
            >::overflow (
                typename basic_bz2_ostreambuf<
                    Elem,Tr,ElemA,ByteT,ByteAT
                    >::int_type c
                )
{
    int w = static_cast<int>(this->pptr() - this->pbase());
    if (c != EOF) {
            *this->pptr() = c;
            ++w;
        }
        if ( bzip2_to_stream( this->pbase(), w)) {
            this->setp( this->pbase(), this->epptr());
            return c;
        } else
            return EOF;
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
bool basic_bz2_ostreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::bzip2_to_stream(
        typename basic_bz2_ostreambuf<
            Elem,Tr,ElemA,ByteT,ByteAT
            >::char_type* buffer_,
        std::streamsize buffer_size_
        )
{
    std::streamsize written_byte_size=0, total_written_byte_size = 0;

    m_bzip2_stream.next_in=(byte_buffer_type)buffer_;
    m_bzip2_stream.avail_in=buffer_size_*sizeof(char_type);
    m_bzip2_stream.avail_out=static_cast<unsigned int>(m_output_buffer.size());
    m_bzip2_stream.next_out=&(m_output_buffer[0]);
    size_t remainder=0;

    do
    {
        m_err = BZ2_bzCompress (&m_bzip2_stream, BZ_RUN );

        if (m_err == BZ_RUN_OK  || m_err == BZ_STREAM_END)
        {
            written_byte_size= static_cast<std::streamsize>(m_output_buffer.size()) - m_bzip2_stream.avail_out;
            total_written_byte_size+=written_byte_size;
            // output buffer is full, dumping to ostream
            m_ostream.write(
                (const char_type*) &(m_output_buffer[0]),
                static_cast<std::streamsize>( written_byte_size/sizeof(char_type) )
                );

            // checking if some bytes were not written.
            if ( (remainder = written_byte_size%sizeof(char_type))!=0)
            {
                // copy to the beginning of the stream
                std::memmove(
                    &(m_output_buffer[0]),
                    &(m_output_buffer[written_byte_size-remainder]),
                    remainder);

            }

            m_bzip2_stream.avail_out=static_cast<unsigned int>(m_output_buffer.size()-remainder);
            m_bzip2_stream.next_out=&m_output_buffer[remainder];
        }
    }
    while (m_bzip2_stream.avail_in != 0 && m_err == BZ_RUN_OK);

    return m_err == BZ_RUN_OK || m_err == BZ_FLUSH_OK;
}

template<
    typename Elem,
    typename Tr,
    typename ElemA,
    typename ByteT,
    typename ByteAT
>
std::streamsize basic_bz2_ostreambuf<
    Elem,Tr,ElemA,ByteT,ByteAT
    >::flush(int flush_mode)
{
    std::streamsize written_byte_size=0, total_written_byte_size=0;

    int const buffer_size = static_cast< int >( pptr() - pbase() ); // amount of data currently in buffer

    m_bzip2_stream.next_in=(byte_buffer_type)pbase();
    m_bzip2_stream.avail_in=static_cast< unsigned int >(buffer_size*sizeof(char_type));
    m_bzip2_stream.avail_out=static_cast< unsigned int >(m_output_buffer.size());
    m_bzip2_stream.next_out=&(m_output_buffer[0]);
    size_t remainder=0;

    do
    {
        m_err = BZ2_bzCompress (&m_bzip2_stream, flush_mode);
        if (m_err == BZ_FINISH_OK || m_err == BZ_STREAM_END)
        {
            written_byte_size=
                static_cast<std::streamsize>(m_output_buffer.size())
                - m_bzip2_stream.avail_out;
            total_written_byte_size+=written_byte_size;
            // output buffer is full, dumping to ostream
            m_ostream.write(
                (const char_type*) &(m_output_buffer[0]),
                static_cast<std::streamsize>( written_byte_size/sizeof(char_type)*sizeof(char) )
                );

            // checking if some bytes were not written.
            if ( (remainder = written_byte_size%sizeof(char_type))!=0)
            {
                // copy to the beginning of the stream
                std::memmove(
                    &(m_output_buffer[0]),
                    &(m_output_buffer[written_byte_size-remainder]),
                    remainder);

            }

            m_bzip2_stream.avail_out=static_cast<unsigned int>(m_output_buffer.size()-remainder);
            m_bzip2_stream.next_out=&(m_output_buffer[remainder]);
        }
    } while (m_err == BZ_FINISH_OK);

    m_ostream.flush();

    return total_written_byte_size;
}

// --------------------------------------------------------------------------
// Class basic_bz2_ostreambase
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bz2_ostreambase : virtual public std::basic_ios<Elem,Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr>& ostream_reference;
    typedef basic_bz2_ostreambuf<
        Elem,Tr,ElemA,ByteT,ByteAT> bzip2_streambuf_type;

    basic_bz2_ostreambase(
        ostream_reference ostream_,
        size_t block_size_100k_ ,
        size_t verbosity_ ,
        size_t work_factor_,
        size_t buffer_size_
        )
        : m_buf(ostream_,block_size_100k_, verbosity_, work_factor_, buffer_size_)
    {
        this->init(&m_buf );
    };

    bzip2_streambuf_type* rdbuf() { return &m_buf; };

private:
    bzip2_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_bz2_ostream
// --------------------------------------------------------------------------

template<
    typename Elem,
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
>
class basic_bz2_ostream :
    public basic_bz2_ostreambase<Elem,Tr,ElemA,ByteT,ByteAT>,
    public std::basic_ostream<Elem,Tr>
{
public:
    typedef basic_bz2_ostreambase<
        Elem,Tr,ElemA,ByteT,ByteAT> bzip2_ostreambase_type;
    typedef std::basic_ostream<Elem,Tr> ostream_type;
    typedef ostream_type& ostream_reference;

    basic_bz2_ostream(
        ostream_reference ostream_,
        size_t block_size_100k_ = 9,
        size_t verbosity_ = 0,
        size_t work_factor_ = 30,
        size_t buffer_size_ = BZ2_OUTPUT_DEFAULT_BUFFER_SIZE
        )
    :
        bzip2_ostreambase_type(ostream_,block_size_100k_, verbosity_, work_factor_,buffer_size_),
        ostream_type(bzip2_ostreambase_type::rdbuf())
    {

    };

    basic_bz2_ostream& add_header();
    basic_bz2_ostream& zflush()
    {
        this->flush(); this->rdbuf()->flush(); return *this;
    };

#ifdef _WIN32
private:
    void _Add_vtordisp1() { } // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() { } // Required to avoid VC++ warning C4250
#endif
};

// --------------------------------------------------------------------------
// Typedefs
// --------------------------------------------------------------------------

typedef basic_bz2_ostream<char>    bz2_ostream;
typedef basic_bz2_ostream<wchar_t> bz2_wostream;

} // namespace seqan3::contrib
