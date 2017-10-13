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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#pragma once

#include <cassert>
#include <streambuf>
#include <istream>
#include <ostream>

#include <range/v3/iterator_range.hpp>

#include <seqan3/core/concept/core.hpp>
#include <seqan3/io/detail/chunking.hpp>
#include <seqan3/io/concept/stream.hpp>

namespace seqan3::detail
{

// ============================================================================
// Tags
// ============================================================================

using input_direction = std::true_type;

using output_direction = std::false_type;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class StreamBuffer
// ----------------------------------------------------------------------------

/*!
 * @class StreamBuffer
 * @headerfile <seqan/stream.h>
 * @brief Buffer to use in stream.
 *
 * @signature template <typename TValue[, typenam TTraits]>
 *            class StreamBuffer : public std::basic_streambuf<TValue, TTraits>;
 *
 * @tparam TValue  The value type of the stream buffer.
 * @tparam TTraits The traits to use, defaults to <tt>std::char_traits&lt;TValue&gt;</tt>.
 */

// TODO(holtgrew): Add documentation for member functions.

// Unfortunately some of the most useful members of basic_streambuf are
// protected, so we define a subclass to cast and access them

template <typename char_t,
          typename traits_t = std::char_traits<char_t>>
class stream_buffer : public std::basic_streambuf<char_t, traits_t>
{
    using base_type = std::basic_streambuf<char_t, traits_t>;
public:
    using base_type::eback;
    using base_type::gptr;
    using base_type::egptr;

    using base_type::pbase;
    using base_type::pptr;
    using base_type::epptr;

    template <typename offset_t, typename direction_t>
    void advance_chunk(offset_t const ofs, direction_t const & /*dir*/)
    {
        if constexpr (std::is_same_v<direction_t, input_direction>)
            this->gbump(ofs);
        else
            this->pbump(ofs);
    }

    template <typename direction_t>
    void next_chunk(direction_t const & /*dir*/)
    {
        if constexpr (std::is_same_v<direction_t, input_direction>)
        {
            if (gptr() == egptr())
                this->underflow();
        }
        else
        {
            if (pptr() == epptr())
                this->overflow(base_type::traits_type::eof());
        }
    }

    template <typename offset_t, typename direction_t>
    void advance(offset_t const ofs, direction_t const & dir)
    {
        assert(ofs >= 0);
        advanve_impl(ofs, dir);
    }

protected:
    template <typename direction_t>
    size_t chunk_size(direction_t const & /*dir*/)
    {
        if constexpr (std::is_same_v<direction_t, input_direction>)
            return egptr() - gptr();
        else
            return epptr() - pptr();
    }

    template <typename TOffset>
    typename std::streampos
    seekoff(TOffset ofs, std::ios_base::seekdir way, std::ios_base::openmode which)
    {
        return base_type::seekoff(ofs, way, which);
    }

    template <typename offset_t, typename direction_t>
    void advanve_impl(offset_t const /*ofs*/, direction_t const & /*dir*/);
};

template <typename char_t, typename traits_t>
template <typename offset_t, typename direction_t>
void stream_buffer<char_t, traits_t>::advanve_impl(offset_t ofs, direction_t const & dir)
{
    size_t left = chunk_size(dir);
    if (/*SEQAN_LIKELY*/((size_t)ofs <= left))
    {
        advance_chunk(ofs, dir);
        return;
    }

    while (true)
    {
        size_t adv = std::min(static_cast<size_t>(ofs), left);
        advance_chunk(adv, dir);
        ofs -= adv;
        if (ofs == 0)
            return;

        next_chunk(dir);
        left = chunk_size(dir);

        if (/*SEQAN_UNLIKELY*/(left == 0))
        {
            // if chunking isn't available try to seek
            typename traits_t::pos_type res =
                seekoff(ofs,
                        std::ios_base::cur,
                        (std::is_same_v<direction_t, input_direction>) ? std::ios_base::in : std::ios_base::out);

            // if seek doesn't work manually skip characters (when reading)
            if (res == typename traits_t::pos_type(typename traits_t::off_type(-1)))
            {
                if constexpr (std::is_same_v<direction_t, input_direction>)
                {
                    for (; ofs != 0; --ofs)
                        this->sbumpc();
                }
                else
                {
                    for (; ofs != 0; --ofs)
                        this->sputc('\0');
                }
            }
            return;
        }
    }
}

// ----------------------------------------------------------------------------
// Class StreamIterator
// ----------------------------------------------------------------------------

/*!
 * @class StreamIterator
 * @extends Iter
 * @brief Abstract base class for input and output stream iterators.
 *
 * @signature template <typename TStream, typename TDirection>
 *            class Iter<TStream, StreamIterator<TDirection> >;
 *
 * @tparam TStream    The @link StreamConcept @endlink to iterate over.
 * @tparam TDirection The iterator direction, one of the @link DirectionTags @endlink.
 */

// ----------------------------------------------------------------------------
// Class Input StreamIterator
// ----------------------------------------------------------------------------

/*!
 * @class InputStreamIterator Input StreamIterator
 * @extends StreamIterator
 * @brief @link Iter @endlink specialiazion for reading from @link StreamConcept streams @endlink.
 *
 * @signature template <typename TStream>
 *            class Iter<TStream, StreamIterator<Input> >;
 *
 * @tparam TStream    The @link StreamConcept @endlink to iterate over.
 */

template <typename stream_t>
//    requires input_stream_concept<stream_t>
class istream_chunk_adaptor_iterator : public chunk_decorator<istream_chunk_adaptor_iterator<stream_t>>
{
    friend class chunk_decorator<istream_chunk_adaptor_iterator<stream_t>>;

    using chunk_base_type    = chunk_decorator<istream_chunk_adaptor_iterator<stream_t>>;
    using basicbuf_type      = std::basic_streambuf<typename stream_t::char_type, typename stream_t::traits_type>;
public:

    // Forwarding char_type and traits_type
    using char_type          = typename stream_t::char_type;
    using traits_type        = typename stream_t::traits_type;
    using streambuf_type     = stream_buffer<char_type, traits_type>;
    using istream_type       = std::basic_istream<char_type, traits_type>;

    // Global member types.
    using value_type         = char_type;
    using difference_type    = typename traits_type::off_type;
    using pointer            = std::add_pointer_t<value_type>;
    using reference          = char_type;
    using iterator_category  = std::input_iterator_tag;

    /*!
     * @fn InputStreamIterator::Iter
     * @brief The constructors.
     *
     * @signature Iter::Iter();
     * @signature Iter::Iter(stream);
     * @signature Iter::Iter(streamBuffer);
     *
     * @param[in] stream    The <tt>TStream</tt> to read from.
     * @param[in] streamBuf A @link StreamBuffer @endlink to read from.
     *
     * Allows default construction, construction from stream, as well as from a @link StreamBuffer @endlink.
     */
    constexpr istream_chunk_adaptor_iterator() = default;

    // Requirements: the input stream_type must fulfil some expectations.
    istream_chunk_adaptor_iterator(istream_type & stream) :
                           streambuf_ptr(static_cast<streambuf_type *>(stream.rdbuf()))
    {
        stream.exceptions(std::ios_base::badbit);
    }

    istream_chunk_adaptor_iterator(basicbuf_type * buf) :
                           streambuf_ptr(static_cast<streambuf_type *>(buf))
    {}

    istream_chunk_adaptor_iterator(istream_chunk_adaptor_iterator const & /*other*/) = default;
    istream_chunk_adaptor_iterator(istream_chunk_adaptor_iterator && /*other*/) = default;

    istream_chunk_adaptor_iterator & operator=(istream_chunk_adaptor_iterator const & /*other*/) = default;
    istream_chunk_adaptor_iterator & operator=(istream_chunk_adaptor_iterator && /*other*/) = default;

    ~istream_chunk_adaptor_iterator() = default;

    // operator implementation.
    reference
    operator*() const
    {
        return streambuf_ptr->sgetc();
    }

    istream_chunk_adaptor_iterator &
    operator++(/*pre*/)
    {
        streambuf_ptr->sbumpc();
        return *this;
    }

    istream_chunk_adaptor_iterator
    operator++(int /*post*/)
    {
        auto tmp{*this};
        ++(*this);
        return tmp;
    }

    inline bool
    operator==(istream_chunk_adaptor_iterator const & rhs) const
    {
        return equal(rhs);
    }

    inline bool
    operator!=(istream_chunk_adaptor_iterator const & rhs) const
    {
        return !equal(rhs);
    }

    template <typename intergral_t>
    inline void
    advance_stream(intergral_t const offset)
    {
        streambuf_ptr->advance(offset, input_direction{});
    }

private:
    streambuf_type * streambuf_ptr{nullptr};

    bool
    equal(istream_chunk_adaptor_iterator const & other) const
    {
        return at_eof() == other.at_eof();
    }

    bool
    at_eof() const
    {
        if (/*SEQAN_UNLIKELY*/(streambuf_ptr))
        {
            if (/*SEQAN_LIKELY*/(streambuf_ptr->gptr() < streambuf_ptr->egptr()))
                return true;
            else
                return !traits_type::eq_int_type(streambuf_ptr->sgetc(), traits_type::eof());
        }
        return false;
    }

    auto chunk_current() const
    {
        return streambuf_ptr->gptr();
    }

    auto chunk_end() const
    {
        return streambuf_ptr->egptr();
    }

    template <typename intergral_t>
    inline void
    next_chunk_impl(intergral_t const /*offset*/)
    {
        streambuf_ptr->next_chunk(input_direction{});
    }

    template <typename offset_t>
    inline void
    advance_chunk_impl(offset_t const offset)
    {
        streambuf_ptr->advance_chunk(offset, input_direction{});
    }

    inline void
    trim_trailing_impl() noexcept
    {
        // no-op
    }
};

// ----------------------------------------------------------------------------
// Class ostream_chunk_adaptor_iterator
// ----------------------------------------------------------------------------

/*!
 * @class OutputStreamIterator Output StreamIterator
 * @extends StreamIterator
 * @brief @link Iter @endlink specialiazion for writing to @link StreamConcept streams @endlink.
 *
 * @signature template <typename TStream>
 *            class Iter<TStream, StreamIterator<Output> >;
 *
 * @tparam TStream    The @link StreamConcept @endlink to iterate over.
 */
template <typename stream_t>
//    requires onput_stream_concept<stream_t>
class ostream_chunk_adaptor_iterator : public chunk_decorator<ostream_chunk_adaptor_iterator<stream_t>>
{
    friend class chunk_decorator<ostream_chunk_adaptor_iterator<stream_t>>;

    using basicbuf_type      = std::basic_streambuf<typename stream_t::char_type, typename stream_t::traits_type>;
    using chunk_base_type    = chunk_decorator<istream_chunk_adaptor_iterator<stream_t>>;

public:

    // Forwarding char_type and traits_type
    using char_type          = typename stream_t::char_type;
    using traits_type        = typename stream_t::traits_type;
    using streambuf_type     = stream_buffer<char_type, traits_type>;
    using ostream_type       = std::basic_ostream<char_type, traits_type>;

    // Global member types.
    using value_type         = void;
    using difference_type    = std::ptrdiff_t;  // required for the output_iterator_concept
    using pointer            = void;
    using reference          = void;
    using iterator_category  = std::output_iterator_tag;

    /*!
     * @fn Iter::Iter
     * @brief Constructor.
     *
     * @signature Iter::Iter()
     * @signature Iter::Iter(stream)
     * @signature Iter::Iter(streamBuf)
     *
     * @param[in] stream    The <tt>TStream</tt> to write to.
     * @param[in] streamBuf A @link StreamBuffer @endlink to write to.
     *
     * Allows default construction, construction from stream, as well as from a @link StreamBuffer @endlink.
     */
    ostream_chunk_adaptor_iterator() = default;

    ostream_chunk_adaptor_iterator(ostream_type & stream) :
                           streambuf_ptr(static_cast<streambuf_type *>(stream.rdbuf()))
    {
        stream.exceptions(std::ios_base::badbit);
    }

    ostream_chunk_adaptor_iterator(basicbuf_type * buf) :
                           streambuf_ptr(static_cast<streambuf_type *>(buf))
    {}

    ostream_chunk_adaptor_iterator(ostream_chunk_adaptor_iterator const & /*other*/) = default;
    ostream_chunk_adaptor_iterator(ostream_chunk_adaptor_iterator && /*other*/) = default;

    ostream_chunk_adaptor_iterator & operator=(ostream_chunk_adaptor_iterator const & /*other*/) = default;
    ostream_chunk_adaptor_iterator & operator=(ostream_chunk_adaptor_iterator && /*other*/) = default;

    ~ostream_chunk_adaptor_iterator() = default;

    template <typename value_t>
        requires convertible_to_concept<value_t, char_type>
    ostream_chunk_adaptor_iterator & operator=(value_t const & val)
    {
        streambuf_ptr->sputc(val);
        return *this;
    }

    ostream_chunk_adaptor_iterator & operator*()
    {
        return *this;
    }

    ostream_chunk_adaptor_iterator & operator++()
    {
        return *this;
    }

    ostream_chunk_adaptor_iterator & operator++(int /*post*/)
    {
        return *this;
    }

    template <typename intergral_t>
    inline void
    advance_stream(intergral_t const offset)
    {
        streambuf_ptr->advance(offset, output_direction{});
    }

private:

    streambuf_type * streambuf_ptr{nullptr};

    auto chunk_current() const
    {
        return streambuf_ptr->pptr();
    }

    auto chunk_end() const
    {
        return streambuf_ptr->epptr();
    }

    template <typename intergral_t>
    inline void
    next_chunk_impl(intergral_t const /*offset*/)
    {
        streambuf_ptr->next_chunk(output_direction{});
    }

    template <typename offset_t>
    inline void
    advance_chunk_impl(offset_t const offset)
    {
        streambuf_ptr->advance_chunk(offset, output_direction{});
    }

    inline void
    trim_trailing_impl() noexcept
    {
        // no-op
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function input_iterator()
// ----------------------------------------------------------------------------

/*!
 * @fn StreamConcept#directionIterator
 * @brief Returns direction iterator for Stream.
 *
 * @signature TDirIter directionIterator(stream, dirTag);
 *
 * @param[in] stream The @link StreamConcept @endlink object to compute iterator for.
 * @param[in] dirTag Direction tag, one of the @link DirectionTags @endlink.
 */

// Returns input itertor pointing to begin of container.
template <typename stream_t>
    requires input_stream_concept<std::decay_t<stream_t>>
inline auto
make_preferred_input_iterator_range(stream_t & stream)
    requires input_stream_concept<std::decay_t<stream_t>>
{
    return std::tuple{istream_chunk_adaptor_iterator<stream_t>{stream}, istream_chunk_adaptor_iterator<stream_t>{}};
}

// Returns back_insert iterator for containers supporting *.push_back().
template <typename stream_t>
    requires output_stream_concept<std::decay_t<stream_t>>
inline auto
make_preferred_output_iterator(stream_t & stream)
{
    return ostream_chunk_adaptor_iterator<stream_t>{stream};
}

// ----------------------------------------------------------------------------
// Function advance()
// ----------------------------------------------------------------------------

template <typename stream_t, typename offset_t>
inline void
advance(istream_chunk_adaptor_iterator<stream_t> & iter,
        offset_t const ofs)
{
    assert(ofs > 0);
    iter.advance_stream(ofs);
}

template <typename stream_t, typename offset_t>
inline void
advance(ostream_chunk_adaptor_iterator<stream_t> & iter,
        offset_t const ofs)
{
    assert(ofs > 0);
    iter.advance_stream(ofs);
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------
// TODO(rrahn): Add when needed
//template <typename TStream, typename TDirection>
//inline typename Position<Iter<TStream, StreamIterator<TDirection> > const>::Type
//position(Iter<TStream, StreamIterator<TDirection> > const & iter)
//{
//    SEQAN_ASSERT(iter.streamBuf != NULL);
//    return iter.streamBuf->pubseekoff(0, std::ios_base::cur,
//                                      (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
//}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------
// TODO(rrahn): Add when needed
//template <typename TStream, typename TDirection, typename TPosition>
//inline void
//setPosition(Iter<TStream, StreamIterator<TDirection> > const & iter, TPosition pos)
//{
//    SEQAN_ASSERT(iter.streamBuf != NULL);
//    iter.streamBuf->pubseekpos(pos, (IsSameType<TDirection, Input>::VALUE)? std::ios_base::in: std::ios_base::out);
//}

}  // namespace seqan3::detail

#ifndef NDEBUG

#include <sstream>

#include <seqan3/core/concept/iterator.hpp>

namespace seqan3::detail
{

static_assert(stream_concept<std::istream>);
static_assert(input_stream_concept<std::istream>);
static_assert(!output_stream_concept<std::istream>);
static_assert(stream_concept<std::ostream>);
static_assert(!input_stream_concept<std::ostream>);
static_assert(output_stream_concept<std::ostream>);

static_assert(input_iterator_concept<istream_chunk_adaptor_iterator<std::stringstream>>);
static_assert(!output_iterator_concept<istream_chunk_adaptor_iterator<std::stringstream>, char>);
static_assert(!forward_iterator_concept<istream_chunk_adaptor_iterator<std::stringstream>>);
static_assert(!bidirectional_iterator_concept<istream_chunk_adaptor_iterator<std::stringstream>>);
static_assert(!random_access_iterator_concept<istream_chunk_adaptor_iterator<std::stringstream>>);

static_assert(!input_iterator_concept<ostream_chunk_adaptor_iterator<std::stringstream>>);
static_assert(output_iterator_concept<ostream_chunk_adaptor_iterator<std::stringstream>, typename std::stringstream::char_type>);
static_assert(!forward_iterator_concept<ostream_chunk_adaptor_iterator<std::stringstream>>);
static_assert(!bidirectional_iterator_concept<ostream_chunk_adaptor_iterator<std::stringstream>>);
static_assert(!random_access_iterator_concept<ostream_chunk_adaptor_iterator<std::stringstream>>);

} // namespace seqan3::detail

#endif // NDEBUG
