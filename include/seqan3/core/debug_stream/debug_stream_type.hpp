// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <iosfwd>
#include <stdexcept>

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/debug_stream/default_printer.hpp>

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::fmtflags2
// ------------------------------------------------------------------

//!\brief Flags that change the behaviour of the seqan3::debug_stream.
//!\ingroup core_debug_stream
//!\implements seqan3::enum_bitwise_operators
//!\sa seqan3::enum_bitwise_operators enables combining enum values.
enum fmtflags2
{
    none = 0,                     //!< No flag is set.
    utf8 = 1,                     //!< Enables use of non-ASCII UTF8 characters in formatted output.
    small_int_as_number = 1 << 1, //!< `int8_t` and `uint8_t` are often aliases for `signed char` and `unsigned char`
                                  //!< resp. resulting in chars being printed; this options prints them as numbers.
    default_ = small_int_as_number
};

//!\cond DEV
//!\brief Overload bitwise operators for seqan3::fmtflags2.
//!\sa seqan3::enum_bitwise_operators enables combining enum values.
template <>
inline constexpr bool add_enum_bitwise_operators<fmtflags2> = true;
//!\endcond

// ------------------------------------------------------------------
// seqan3::debug_stream_type
// ------------------------------------------------------------------

/*!\brief A "pretty printer" for most SeqAn data structures and related types.
 * \ingroup core_debug_stream
 * \details
 *
 * A global instance of this type exists as seqan3::debug_stream. You can stream to it as you would to std::cout or
 * std::cerr, but the debug stream has special
 * overloads that make certain types streamable (that are not streamable to std::cout). Additionally some
 * data structures are visualised more elaborately via the debug stream and there are extra flags to configure it
 * (seqan3::fmtflags2).
 *
 * \see core_debug_stream
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
     * \brief The standard functions are explicitly set to default.
     * \{
     */
    debug_stream_type() = default;                                      //!< Defaulted
    debug_stream_type(debug_stream_type const &) = default;             //!< Defaulted
    debug_stream_type(debug_stream_type &&) = default;                  //!< Defaulted
    debug_stream_type & operator=(debug_stream_type const &) = default; //!< Defaulted
    debug_stream_type & operator=(debug_stream_type &&) = default;      //!< Defaulted
    ~debug_stream_type() = default;                                     //!< Defaulted

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
    template <typename other_char_t, typename t>
    friend debug_stream_type<other_char_t> & operator<<(debug_stream_type<other_char_t> & s, t && v)
    {
        using t_ = std::remove_cvref_t<t>;

        if constexpr (default_printer::is_printable<t_>)
        {
            default_printer::print(s, v);
        }
        else
        {
            std::string const msg = std::string{"debug_stream has no print overload for type: "} + typeid(v).name();
            throw std::runtime_error{msg};
        }
        return s;
    }

    //!\brief This overloads enables forwarding std::endl and other manipulators.
    debug_stream_type & operator<<(std::ostream & (*fp)(std::ostream &))
    {
        *stream << fp;
        return *this;
    }

    //!\}

    template <typename T>
    friend struct debug_stream_printer;

    template <typename T>
    friend struct std_printer;

    //!\brief This type is std::ios_base::fmtflags
    using fmtflags = typename std::basic_ostream<char_t>::fmtflags;

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

    //!\}

private:
    //!\brief Pointer to the output stream.
    std::basic_ostream<char_t> * stream /* = &std::cerr*/;

    //!\brief The SeqAn specific flags to the stream.
    fmtflags2 flgs2{fmtflags2::default_};
};

template <typename value_t>
    requires (std::is_same_v<std::remove_cvref_t<value_t>, int8_t>
              || std::is_same_v<std::remove_cvref_t<value_t>, uint8_t>
              || std::is_same_v<std::remove_cvref_t<value_t>, fmtflags2>)
struct debug_stream_printer<value_t>
{
    struct print_fn
    {
        template <typename stream_t>
        constexpr void operator()(stream_t & s, int8_t const v) const
        {
            if ((s.flags2() & fmtflags2::small_int_as_number) == fmtflags2::small_int_as_number)
                *s.stream << static_cast<int>(v);
            else
                *s.stream << v;
        }

        template <typename stream_t>
        constexpr void operator()(stream_t & s, uint8_t const v) const
        {
            if ((s.flags2() & fmtflags2::small_int_as_number) == fmtflags2::small_int_as_number)
                *s.stream << static_cast<unsigned>(v);
            else
                *s.stream << v;
        }

        template <typename stream_t>
        constexpr void operator()(stream_t & s, fmtflags2 const flag) const
        {
            s.setf(flag);
        }
    };

    static constexpr print_fn print{};
};

} // namespace seqan3
