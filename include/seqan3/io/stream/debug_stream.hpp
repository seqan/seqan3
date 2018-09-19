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

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <ios>
#include <ostream>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::fmtflags2
// ------------------------------------------------------------------

//!\brief Flags that change the behaviour of the seqan3::debug_stream.
//!\ingroup stream
enum fmtflags2
{
    none                = 0,        //!< No flag is set.
    utf8                = 1,        //!< Enables use of non-ASCII UTF8 characters in formatted output.
    small_int_as_number = 1 << 1,   //!< `int8_t` and `uint8_t` are often aliases for `signed char` and `unsigned char`
                                    //!< resp. resulting in chars being printed; this options prints them as numbers.
    default_            = small_int_as_number
};

//!\brief Overload bitwise operators for seqan3::fmtflags2.
template <>
constexpr bool add_enum_bitwise_operators<fmtflags2> = true;

// ------------------------------------------------------------------
// seqan3::debug_stream_type
// ------------------------------------------------------------------

/*!\brief A "pretty printer" for most SeqAn data structures and related types.
 * \ingroup stream
 * \details
 *
 * A global instance of this type exists as seqan3::debug_stream. You can stream to it as you would to std::cout or
 * std::cerr, but the debug stream has special
 * overloads that make certain types streamable (that are not streamable to std::cout). Additionally some
 * data structures are visualised more elaborately via the debug stream and there are extra flags to configure it
 * (seqan3::fmtflags2).
 *
 * ## Example
 *
 *  Simple usage:
 *
 * \snippet test/snippet/io/stream/debug_stream.cpp usage
 *
 * Changing flags:
 *
 * \snippet test/snippet/io/stream/debug_stream.cpp flags
 *
 * See seqan3::fmtflags2 for more details.
 *
 * \attention This class does not yet model seqan3::ostream_concept fully, \todo implement.
 */
class debug_stream_type
{
public:
    /*!\name Constructor, destructor and assignment.
     * \brief The standard functions are explicitly defaulted.
     * \{
     */
    constexpr debug_stream_type() = default;
    constexpr debug_stream_type(debug_stream_type const &) = default;
    constexpr debug_stream_type(debug_stream_type &&) = default;
    constexpr debug_stream_type & operator= (debug_stream_type const &) = default;
    constexpr debug_stream_type & operator= (debug_stream_type &&) = default;
    ~debug_stream_type() = default;

    //!\brief Construction from an output stream.
    constexpr explicit debug_stream_type(std::ostream & out) : stream{&out}
    {}
    //!\}

    /*!\name Miscelleneous
     * \{
     */
    /*!\brief Change the underlying output stream.
     * \param out Reference to the new output stream.
     *
     * \details
     *
     * The actual underlying stream that is printed to defaults to std::cerr, but can be changed via this function.
     * You can set any kind of output stream, e.g. a std::ostringstream or a
     * std::ofstream if you want to write to a file, but please be aware that the debug_stream never takes ownership of
     * the underlying stream so you need to take special care that its object lifetime does not end before the
     * debug_stream's.
     *
     * \snippet test/snippet/io/stream/debug_stream.cpp set_underlying_stream
     *
     * In the case where you wish to print to some stream object locally, instead create you own debug stream:
     *
     * \snippet test/snippet/io/stream/debug_stream.cpp set_underlying_stream2
     */
    void set_underlying_stream(std::ostream & out)
    {
        stream = &out;
    }
    //!\}

    /*!\name Formatted output
     * \{
     */
    //!\brief Forwards to the underlying stream object.
    template <typename t>
    debug_stream_type & operator<<(t const & v)
    {
        *stream << v;
        return *this;
    }

    //!\cond
    debug_stream_type & operator<<(int8_t const v)
    {
        if ((flags2() & fmtflags2::small_int_as_number) == fmtflags2::small_int_as_number)
            *stream << static_cast<int>(v);
        else
            *stream << v;
        return *this;
    }

    debug_stream_type & operator<<(uint8_t const v)
    {
        if ((flags2() & fmtflags2::small_int_as_number) == fmtflags2::small_int_as_number)
            *stream << static_cast<unsigned>(v);
        else
            *stream << v;
        return *this;
    }
    //!\endcond
    //!\}

    /*!\name Format flags (std::ios_base::fmtflags)
     * \brief std::ios_base::fmtflags that modify the stream's behaviour.
     * \{
     */
    //!\brief Retrieve the format flags from the stream.
    std::ios_base::fmtflags flags() const
    {
        return stream->flags();
    }

    //!\brief Replace the current flags on the stream with the given argument.
    std::ios_base::fmtflags flags(std::ios_base::fmtflags const flgs)
    {
        return stream->flags(flgs);
    }

    //!\brief Set the format flag(s) on the stream (current flags are ORed with the argument).
    void setf(std::ios_base::fmtflags const flag)
    {
        stream->setf(flag);
    }

    //!\brief Unset the format flag(s) on the stream.
    void unsetf(std::ios_base::fmtflags const flag)
    {
        stream->unsetf(flag);
    }

    //!\copybrief setf()
    debug_stream_type & operator<<(std::ios_base::fmtflags const flag)
    {
        setf(flag);
        return *this;
    }
    //!\}

    /*!\name Format flags  (seqan3::fmtflags2)
     * \brief SeqAn specific debug flags for the debug stream.
     * \{
     */
    //!\copybrief flags()
    fmtflags2 flags2() const
    {
        return flgs2;
    }

    //!\copybrief flags(std::ios_base::fmtflags const flgs)
    fmtflags2 flags2(fmtflags2 flgs)
    {
        flgs2 = flgs;
        return flgs2;
    }

    //!\copybrief setf()
    void setf(fmtflags2 const flag)
    {
        flgs2 |= flag;
    }

    //!\copybrief unsetf()
    void unsetf(fmtflags2 const flag)
    {
        flgs2 &= ~flag;
    }

    //!\copybrief setf()
    debug_stream_type & operator<<(fmtflags2 const flag)
    {
        setf(flag);
        return *this;
    }
    //!\}

private:
    //!\brief Pointer to the output stream.
    std::ostream *stream = &std::cerr;

    //!\brief The SeqAn specific flags to the stream.
    fmtflags2 flgs2{fmtflags2::default_};
};

// ------------------------------------------------------------------
// seqan3::debug_stream
// ------------------------------------------------------------------

//!\brief A global instance of seqan3::debug_stream_type.
//!\ingroup stream
inline debug_stream_type debug_stream{};

// ------------------------------------------------------------------
// overloads
// ------------------------------------------------------------------

/*!\name Formatted output overloads
 * \{
 */
/*!\brief All alphabets can be printed to the seqan3::debug_stream by their char representation.
 * \tparam alphabet_t Type of the alphabet to be printed; must model seqan3::alphabet_concept.
 * \param s The seqan3::debug_stream.
 * \param l The alphabet letter.
 * \relates seqan3::debug_stream_type
 */
template <alphabet_concept alphabet_t>
inline debug_stream_type & operator<<(debug_stream_type & s, alphabet_t const l)
//!\cond
    requires !ostream_concept<std::ostream, alphabet_t>
//!\endcond
{
    return s << to_char(l);
}

/*!\brief All input ranges can be printed to the seqan3::debug_stream element-wise (if their elements are printable).
 * \tparam rng_t Type of the range to be printed; must model std::ranges::InputRange.
 * \param s The seqan3::debug_stream.
 * \param r The input range.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * If the element type models seqan3::alphabet_concept (and is not an unsigned integer), the range is printed
 * just as if it were a string, i.e. `std::vector<dna4>{dna4::C, dna4:G, dna4::A}` is printed as "CGA".
 *
 * In all other cases the elements are comma separated and the range is enclosed in brackets, i.e.
 * `std::vector<int>{3, 1, 33, 7}` is printed as "[3,1,33,7]".
 */
template <std::ranges::InputRange rng_t>
inline debug_stream_type & operator<<(debug_stream_type & s, rng_t && r)
//!\cond
    requires requires (reference_t<rng_t> l) { { debug_stream << l }; }
//!\endcond
{
    if constexpr (alphabet_concept<reference_t<rng_t>> &&
                  !detail::is_uint_adaptation_v<remove_cvref_t<reference_t<rng_t>>>)
    {
        for (auto && l : r)
            s << l;
    }
    else
    {
        s << '[';
        auto b = ranges::begin(r);
        auto e = ranges::end(r);
        if (b != e)
        {
            s << *b;
            ++b;
        }
        while (b != e)
        {
            s << ',';
            s << *b;
            ++b;
        }
        s << ']';
    }

    return s;
}

//!\}

} // namespace seqan3
