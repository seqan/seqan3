// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <ios>
#include <ostream>
#include <variant>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/range/shortcuts.hpp>
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
 * # Example
 *
 *  Simple usage:
 *
 * \snippet test/snippet/core/debug_stream.cpp usage
 *
 * Changing flags:
 *
 * \snippet test/snippet/core/debug_stream.cpp flags
 *
 * See seqan3::fmtflags2 for more details.
 *
 * \attention This class does not yet model seqan3::OStream fully, \todo implement.
 */
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
     * \snippet test/snippet/core/debug_stream.cpp set_underlying_stream
     *
     * In the case where you wish to print to some stream object locally, instead create you own debug stream:
     *
     * \snippet test/snippet/core/debug_stream.cpp set_underlying_stream2
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
 * \tparam alphabet_t Type of the alphabet to be printed; must model seqan3::Alphabet.
 * \param s The seqan3::debug_stream.
 * \param l The alphabet letter.
 * \relates seqan3::debug_stream_type
 */
template <Alphabet alphabet_t>
inline debug_stream_type & operator<<(debug_stream_type & s, alphabet_t const l)
//!\cond
    requires !OStream<std::ostream, alphabet_t>
//!\endcond
{
    return s << to_char(l);
}

}

namespace seqan3::detail
{

//!\brief Helper function to print elements of a tuple separately.
template<typename tuple_t, std::size_t ...I>
void print_tuple(debug_stream_type & s, tuple_t && t, std::index_sequence<I...> const &)
{
    s << '(';
    ((s << (I == 0 ? "" : ",") << std::get<I>(t)), ...);
    s << ')';
}

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief All tuples can be printed by printing their elements separately.
 * \tparam tuple_t Type of the tuple to be printed; must model seqan3::TupleLike.
 * \param s The seqan3::debug_stream.
 * \param t The tuple.
 * \relates seqan3::debug_stream_type
 */
template <typename tuple_t>
//!\cond
    requires !std::ranges::InputRange<tuple_t> &&
             !Alphabet<tuple_t> && // exclude alphabet_tuple_base
             TupleLike<remove_cvref_t<tuple_t>>
//!\endcond
inline debug_stream_type & operator<<(debug_stream_type & s, tuple_t && t)
{
    detail::print_tuple(s, std::forward<tuple_t>(t),
                        std::make_index_sequence<std::tuple_size_v<remove_cvref_t<tuple_t>>>{});
    return s;
}

/*!\brief A std::variant can be printed by visiting the stream operator for the corresponding type.
 * \tparam    variant_type The underlying type of the variant.
 * \param[in] s            The seqan3::debug_stream.
 * \param[in] v            The variant.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * Note that in case the variant is valueless(_by_exception), nothing is printed.
 */
template <typename variant_type>
//!\cond
    requires detail::is_type_specialisation_of_v<remove_cvref_t<variant_type>, std::variant>
//!\endcond
inline debug_stream_type & operator<<(debug_stream_type & s, variant_type && v)
{
    if (!v.valueless_by_exception())
        std::visit([&s] (auto && arg) {s << arg;}, v);
    else
        s << "<VALUELESS_VARIANT>";
    return s;
}

/*!\brief A std::optional can be printed by printing its value or nothing if valueless.
 * \tparam    optional_type The underlying type of the optional.
 * \param[in] s             The seqan3::debug_stream.
 * \param[in] arg           The std::optional.
 * \relates seqan3::debug_stream_type
 */
template <typename optional_type>
//!\cond
    requires detail::is_type_specialisation_of_v<remove_cvref_t<optional_type>, std::optional>
//!\endcond
inline debug_stream_type & operator<<(debug_stream_type & s, optional_type && arg)
{
    if (arg.has_value())
        s << *arg;
    else
        s << "<VALUELESS_OPTIONAL>";
    return s;
}

/*!\brief All input ranges can be printed to the seqan3::debug_stream element-wise (if their elements are printable).
 * \tparam rng_t Type of the range to be printed; must model std::ranges::InputRange.
 * \param s The seqan3::debug_stream.
 * \param r The input range.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * If the element type models seqan3::Alphabet (and is not an unsigned integer), the range is printed
 * just as if it were a string, i.e. <tt>std::vector<dna4>{'C'_dna4, 'G'_dna4, 'A'_dna4}</tt> is printed as "CGA".
 *
 * In all other cases the elements are comma separated and the range is enclosed in brackets, i.e.
 * `std::vector<int>{3, 1, 33, 7}` is printed as "[3,1,33,7]".
 */
template <std::ranges::InputRange rng_t>
inline debug_stream_type & operator<<(debug_stream_type & s, rng_t && r)
//!\cond
    requires !std::Same<remove_cvref_t<reference_t<rng_t>>, remove_cvref_t<rng_t>> && // prevent recursive instantiation
             requires (reference_t<rng_t> l) { { debug_stream << l }; } &&
             // exclude null-terminated strings:
             !(std::is_pointer_v<std::decay_t<rng_t>> &&
               std::Same<remove_cvref_t<reference_t<rng_t>>, char>)
//!\endcond
{
    if constexpr (Alphabet<reference_t<rng_t>> &&
                  !detail::is_uint_adaptation_v<remove_cvref_t<reference_t<rng_t>>>)
    {
        for (auto && l : r)
            s << l;
    }
    else
    {
        s << '[';
        auto b = begin(r);
        auto e = end(r);
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
