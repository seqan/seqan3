// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <iosfwd>

#include <seqan3/core/add_enum_bitwise_operators.hpp>

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
 * \ingroup core
 * \details
 *
 * A global instance of this type exists as seqan3::debug_stream. You can stream to it as you would to std::cout or
 * std::cerr, but the debug stream has special
 * overloads that make certain types streamable (that are not streamable to std::cout). Additionally some
 * data structures are visualised more elaborately via the debug stream and there are extra flags to configure it
 * (seqan3::fmtflags2).
 *
 * # Example
 *
 *  Simple usage:
 *
 * \include test/snippet/core/debug_stream_usage.cpp
 *
 * Changing flags:
 *
 * \include test/snippet/core/debug_stream_flags.cpp
 *
 * See seqan3::fmtflags2 for more details.
 *
 * \attention This class does not yet model seqan3::output_stream_over fully, \todo implement.
 */
template <typename char_t = char>
class debug_stream_type
{
public:
    /*!\name Constructor, destructor and assignment.
     * \brief The standard functions are explicitly defaulted.
     * \{
     */
    debug_stream_type() = default;                                       //!< Defaulted
    debug_stream_type(debug_stream_type const &) = default;              //!< Defaulted
    debug_stream_type(debug_stream_type &&) = default;                   //!< Defaulted
    debug_stream_type & operator= (debug_stream_type const &) = default; //!< Defaulted
    debug_stream_type & operator= (debug_stream_type &&) = default;      //!< Defaulted
    ~debug_stream_type() = default;                                      //!< Defaulted

    //!\brief Construction from an output stream.
    constexpr explicit debug_stream_type(std::basic_ostream<char_t> & out) : stream{&out}
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
     * \include test/snippet/core/debug_stream_set_underlying_stream.cpp
     *
     * In the case where you wish to print to some stream object locally, instead create you own debug stream:
     *
     * \include test/snippet/core/debug_stream_set_underlying_stream2.cpp
     */
    void set_underlying_stream(std::basic_ostream<char_t> & out)
    {
        stream = &out;
    }
    //!\}

    /*!\name Formatted output
     * \{
     */
    //!\brief Forwards to the underlying stream object.
    template <typename t>
    debug_stream_type & operator<<(t && v)
    {
        *stream << v;
        return *this;
    }

    //!\brief This overloads enables forwarding std::endl and other manipulators.
    debug_stream_type & operator<<(std::ostream&(*fp)(std::ostream&))
    {
        *stream << fp;
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

    //!\brief This type is std::ios_base::fmtflags
    using fmtflags = decltype(std::declval<std::basic_ostream<char_t>>().flags());

    /*!\name Format flags (std::ios_base::fmtflags)
     * \brief std::ios_base::fmtflags that modify the stream's behaviour.
     * \{
     */
    //!\brief Retrieve the format flags from the stream.
    fmtflags flags() const
    {
        return stream->flags();
    }

    //!\brief Replace the current flags on the stream with the given argument.
    fmtflags flags(fmtflags const flgs)
    {
        return stream->flags(flgs);
    }

    //!\brief Set the format flag(s) on the stream (current flags are ORed with the argument).
    void setf(fmtflags const flag)
    {
        stream->setf(flag);
    }

    //!\brief Unset the format flag(s) on the stream.
    void unsetf(fmtflags const flag)
    {
        stream->unsetf(flag);
    }

    //!\copybrief setf()
    debug_stream_type & operator<<(fmtflags const flag)
    {
        setf(flag);
        return *this;
    }
    //!\}

    /*!\name Format flags (seqan3::fmtflags2)
     * \brief SeqAn specific debug flags for the debug stream.
     * \{
     */
    //!\copybrief flags()
    fmtflags2 flags2() const
    {
        return flgs2;
    }

    //!\copybrief flags(fmtflags const flgs)
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
    std::basic_ostream<char_t> *stream/* = &std::cerr*/;

    //!\brief The SeqAn specific flags to the stream.
    fmtflags2 flgs2{fmtflags2::default_};
};

} // namespace seqan3
