// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Rene Rahn <rene.rahn AT fu-berlin.de>
// ==========================================================================

/*
    Requirements
    work on forward_iterator concept
    direction_iterator over the stream_iterator
    stream_iterator can redirect to underlying buffer
    chunked i/o

    read_until(fwd_iter, ignore_f, stop_f)
*/

#pragma once

#include <cctype>
#include <cstring>
#include <string>
#include <tuple>
#include <stdexcept>
#include <algorithm>

#include <seqan3/core/concept/core.hpp>
#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/core/concept/range.hpp>
#include <seqan3/core/concept/stl_container.hpp>
#include <seqan3/core/meta/associated_types.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/io/detail/direction_iterator.hpp>
#include <seqan3/io/detail/stream_iterator.hpp>

namespace seqan3::detail
{

template <typename t>
struct Value
{
    using Type = typename t::value_type;
};

// ----------------------------------------------------------------------------
// Functor AssertFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor, typename TException, typename TContext = decltype(std::ignore), bool RETURN_VALUE = false>
struct AssertFunctor
{
    TFunctor func;

    AssertFunctor() {}

    AssertFunctor(TFunctor & func) :
    func(func)
    {}

    std::string escapeChar(unsigned char val)
    {
        if (val <= '\r')
        {
            static const char * const escapeCodes[14] = {
                "\\0",  "\\1",  "\\2",  "\\3",  "\\4",  "\\5",  "\\6",  "\\a",
                "\\b",  "\\t",  "\\n",  "\\v",  "\\f",  "\\r" };
            return std::string(escapeCodes[val]);
        }
        else if (' ' <= val && val < 128u)
            return std::string() + (char)val;
        else
        {
            char buffer[6]; // 5 + 1, e.g. "\0xff" + trailing zero
            snprintf(buffer, 6, "\\%#2x", (unsigned)val);
            return std::string(buffer);
        }
    }

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        if (/*SEQAN_UNLIKELY(*/!func(val))/*)*/
            throw TException(std::string("Unexpected character '") + escapeChar(val) + "' found. " +
                             getExceptionMessage(func, TContext()));
        return RETURN_VALUE;
    }
};

// ----------------------------------------------------------------------------
// Functor OrFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor1, typename TFunctor2>
struct OrFunctor
{
    TFunctor1 func1;
    TFunctor2 func2;

    OrFunctor()
    {}

    OrFunctor(TFunctor1 const &func1, TFunctor2 const &func2) :
        func1(func1), func2(func2)
    {}

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        return func1(val) || func2(val);
    }

    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        return func1(val) || func2(val);
    }
};

// ----------------------------------------------------------------------------
// Metafunction ExceptionMessage
// ----------------------------------------------------------------------------

template <typename T, typename TSpec = void>
struct ExceptionMessage
{
    static const std::string VALUE;
};

template <typename T, typename TSpec>
const std::string ExceptionMessage<T, TSpec>::VALUE;

// ----------------------------------------------------------------------------
// Function getExceptionMessage()
// ----------------------------------------------------------------------------

template <typename TFunctor, typename TContext>
inline std::string const &
getExceptionMessage(TFunctor const &, TContext const &)
{
    return ExceptionMessage<TFunctor, TContext>::VALUE;
}

// ----------------------------------------------------------------------------
// Exception parse exceptions()
// ----------------------------------------------------------------------------

using  runtime_error = std::runtime_error;

struct parse_error : runtime_error
{
    template <typename TString>
    parse_error(TString const & message) : runtime_error(message)
    {}
};

struct unexpected_end_error : parse_error
{

    unexpected_end_error() : parse_error("Unexpected end of input.")
    {}
};

// ----------------------------------------------------------------------------
// Exception empty_field_error
// ----------------------------------------------------------------------------

struct empty_field_error : parse_error
{
    empty_field_error(std::string fieldName) : parse_error(fieldName + " field was empty.")
    {}
};

// ============================================================================
// Functors
// ============================================================================

// ----------------------------------------------------------------------------
// Functor IsInAlphabet
// ----------------------------------------------------------------------------

template <typename TValue>
   requires alphabet_concept<TValue>
struct IsInAlphabet
{
    template <typename TInValue>
    bool operator() (TInValue const & inVal) const
    {
        TValue val{};
        from_char(val, inVal);
        return val == std::toupper(inVal);
    }

    template <typename TInValue>
        requires alphabet_concept<TInValue>
    bool operator() (TInValue const & inVal) const
    {
        TValue val{};
        from_integral(val, to_integral(inVal));
        return val == inVal;
    }

    bool operator() (TValue const &) const
    {
        return true;
    }
};

// ----------------------------------------------------------------------------
// Functor IsInRange
// ----------------------------------------------------------------------------

template <char FIRST_CHAR, char LAST_CHAR>
struct IsInRange
{
    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        return FIRST_CHAR <= val && val <= LAST_CHAR;
    }
};

template <char FIRST_CHAR, char LAST_CHAR, typename TContext>
struct ExceptionMessage<IsInRange<FIRST_CHAR, LAST_CHAR>, TContext>
{
    static const std::string VALUE;
};

template <char FIRST_CHAR, char LAST_CHAR, typename TContext>
const std::string ExceptionMessage<IsInRange<FIRST_CHAR, LAST_CHAR>, TContext>::VALUE =
std::string("Character in range'") + FIRST_CHAR + "' to '" + LAST_CHAR + "' expected.";

// ----------------------------------------------------------------------------
// Functor EqualsChar
// ----------------------------------------------------------------------------

template <char VALUE>
struct EqualsChar
{
    template <typename TValue>
    bool operator() (TValue const & val) const
    {
        return val == VALUE;
    }
};

template <char CHAR, typename TContext>
struct ExceptionMessage<EqualsChar<CHAR>, TContext>
{
    static const std::string VALUE;
};

template <char CHAR, typename TContext>
const std::string ExceptionMessage<EqualsChar<CHAR>, TContext>::VALUE = std::string("Character '") + CHAR + "' expected.";

// ----------------------------------------------------------------------------
// Functor EqualsDynamicValue
// ----------------------------------------------------------------------------

template <typename TValue>
struct EqualsDynamicValue
{
    TValue val;

    EqualsDynamicValue(TValue const & val) :
    val(val)
    {}

    template <typename TValue2>
    bool operator() (TValue2 const & v) const
    {
        return v == val;
    }
};

template <typename TValue, typename TContext>
inline std::string const &
getExceptionMessage(EqualsDynamicValue<TValue> const & func, TContext const &)
{
    return std::string("Character '") + func.val + "' expected.";
}

// ----------------------------------------------------------------------------
// Composite Functors
// ----------------------------------------------------------------------------
// Don't use isblank() or isspace() as it they seem to be slower than our functors (due to inlining)

typedef EqualsChar<'\t'>                                        IsTab;
typedef EqualsChar<' '>                                         IsSpace;
typedef OrFunctor<IsSpace, IsTab>                               IsBlank;
typedef OrFunctor<EqualsChar<'\n'>, EqualsChar<'\r'> >          IsNewline;
typedef OrFunctor<IsBlank, IsNewline>                           IsWhitespace;
typedef IsInRange<'!', '~'>                                     IsGraph;
typedef OrFunctor<IsInRange<'a', 'z'>, IsInRange<'A', 'Z'> >    IsAlpha;
typedef IsInRange<'0', '9'>                                     IsDigit;
typedef OrFunctor<IsAlpha, IsDigit>                             IsAlphaNum;

constexpr auto always_true = [](auto &&...){ return true; };
constexpr auto always_false = [](auto &&...){ return false; };

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function put()                                                      [Iter]
// ----------------------------------------------------------------------------

template <typename value_t, typename output_t>
inline void
put(value_t const & val,
    output_t && o_iter)
{
    *o_iter = val;
    ++o_iter;
}

// ----------------------------------------------------------------------------
// Function _write(); Element-wise
// ----------------------------------------------------------------------------

template <typename input_t, typename integral_t, typename output_t>
inline void
write(input_t && i_iter,
      integral_t const n,
      output_t && o_iter)
    requires input_iterator_concept<std::remove_reference_t<input_t>> &&
             output_iterator_concept<std::remove_reference_t<output_t>, value_type_t<std::remove_reference_t<input_t>>> &&
             assignable_concept<value_type_t<std::remove_reference_t<output_t>>&, value_type_t<std::remove_reference_t<input_t>>>
{
    std::cout << "simple copy\n";
    std::copy_n(std::forward<input_t>(i_iter), n, std::forward<output_t>(o_iter));
}

// ----------------------------------------------------------------------------
// Function _write(); Chunked
// ----------------------------------------------------------------------------

// How can we enable chunking for pointers?
// We can get a direction iterator for an pointer type.

// TODO(rrahn): Add later
template <typename input_t, typename integral_t, typename output_t>
inline void
write(input_t && i_iter,
      integral_t n,
      output_t && o_iter)
    requires input_iterator_concept<std::remove_reference_t<input_t>> &&
             std::is_base_of_v<chunk_decorator<std::remove_reference_t<input_t>>, std::remove_reference_t<input_t>> &&
             output_iterator_concept<std::remove_reference_t<output_t>, value_type_t<std::remove_reference_t<input_t>>> &&
             std::is_base_of_v<chunk_decorator<std::remove_reference_t<output_t>>, std::remove_reference_t<output_t>> &&
             assignable_concept<value_type_t<std::remove_reference_t<output_t>>&, value_type_t<std::remove_reference_t<input_t>>>
{
   using target_size_type = size_type_t<output_t>;

   using std::size;
   using std::begin;

   std::cout << "I am here baby!" << std::endl;
   while (n != static_cast<integral_t>(0))
   {
       auto ichunk = i_iter.get_chunk();
       auto ochunk = o_iter.get_chunk();

       target_size_type min_chunk_size = std::min(static_cast<target_size_type>(size(ichunk)),
                                                  static_cast<target_size_type>(size(ochunk)));

       if (/*SEQAN_UNLIKELY*/(min_chunk_size == 0u))
       {
           i_iter.next_chunk(n);
           o_iter.next_chunk(n);
           ichunk = i_iter.get_chunk();
           ochunk = o_iter.get_chunk();
           min_chunk_size = std::min(static_cast<target_size_type>(size(ichunk)),
                                     static_cast<target_size_type>(size(ochunk)));
           if (/*SEQAN_UNLIKELY*/(min_chunk_size == 0u))
           {
               std::copy_n(std::forward<input_t>(i_iter), n, std::forward<output_t>(o_iter));  // fall back to no-chunking version.
               return;
           }
       }

       if (min_chunk_size > static_cast<target_size_type>(n))
           min_chunk_size = static_cast<target_size_type>(n);

       std::copy_n(begin(ichunk), min_chunk_size, begin(ochunk));

       i_iter.advance_chunk(min_chunk_size);
       o_iter.advance_chunk(min_chunk_size);
       n -= min_chunk_size;
   }
}

// chunked, target is pointer (e.g. read_raw_pod)
template <typename input_t, typename integral_t, typename output_t>
inline void
write(input_t && i_iter,
      integral_t n,
      output_t && o_ptr)
    requires input_iterator_concept<std::remove_reference_t<input_t>> &&
             std::is_base_of_v<chunk_decorator<std::remove_reference_t<input_t>>, std::remove_reference_t<input_t>> &&
             std::is_pointer_v<std::remove_reference_t<output_t>> &&
             assignable_concept<value_type_t<std::remove_reference_t<output_t>>&, value_type_t<std::remove_reference_t<input_t>>>
{
    // we need the size type of an iterator?
   using source_size_t = size_type_t<input_t>;

   using std::size;
   using std::begin;

   while (n != static_cast<integral_t>(0))
   {
       auto ichunk = i_iter.get_chunk();
       source_size_t chunk_size = size(ichunk);

       if (/*SEQAN_UNLIKELY*/(chunk_size == 0u))
       {
           next_chunk(i_iter, n);
           ichunk = i_iter.get_chunk();
           source_size_t chunk_size = size(ichunk);
           if (/*SEQAN_UNLIKELY*/(chunk_size == 0u))
           {
               std::copy_n(std::forward<input_t>(i_iter), n, std::forward<output_t>(o_ptr));
               return;
           }
       }

       if (chunk_size > static_cast<source_size_t>(n))
           chunk_size = static_cast<source_size_t>(n);

           std::copy_n(begin(ichunk), chunk_size, std::forward<output_t>(o_ptr));

           i_iter.advance_chunk(chunk_size);                          // advance input iterator
           o_ptr += chunk_size;
           n -= chunk_size;
       }
}

// chunked, source is pointer (e.g. readRawPod)
template <typename input_t, typename integral_t, typename output_t>
inline void
write(input_t && i_ptr,
      integral_t n,
      output_t && o_iter)
    requires std::is_pointer_v<std::remove_reference_t<input_t>> &&
             output_iterator_concept<std::remove_reference_t<output_t>, value_type_t<std::remove_reference_t<input_t>>> &&
             std::is_base_of_v<chunk_decorator<std::remove_reference_t<output_t>>, std::remove_reference_t<output_t>> &&
             assignable_concept<value_type_t<std::remove_reference_t<output_t>>&, value_type_t<std::remove_reference_t<input_t>>>
{
    using output_size_t = size_type_t<output_t>;

    while (n != static_cast<integral_t>(0))
    {
        auto ochunk = get_chunk(o_iter);
        output_size_t chunk_size = size(ochunk);

        if (/*SEQAN_UNLIKELY*/(chunk_size == 0u))
        {
            next_chunk(o_iter, n);
            ochunk = get_chunk(o_iter);
            chunk_size = size(ochunk);
            if (/*SEQAN_UNLIKELY*/(chunk_size == 0u))
            {
                std::copy_n(std::forward<input_t>(i_ptr), n, std::forward<output_t>(o_iter));
               return;
            }
        }

        if (chunk_size > static_cast<output_size_t>(n))
            chunk_size = static_cast<output_size_t>(n);

        std::copy_n(std::forward<input_t>(i_ptr), chunk_size, begin(ochunk));

        i_ptr += chunk_size;                      // advance input iterator
        advance_chunk(o_iter, chunk_size);
        n -= chunk_size;
   }
}

// TODO(rrahn): Do we really need this?
//template <typename TOValue, typename TIValue, typename TSize>
//inline SEQAN_FUNC_ENABLE_IF(And< Is<CharConcept<TOValue> >,
//                            Is<CharConcept<TIValue> > >, void)
//write(TOValue * &optr, TIValue *iptr, TSize n)
//{
//    std::memcpy(optr, iptr, n);
//    optr += n;
//}

// TODO(rrahn): not genric in sense of char and wchar or even bigger char.
// template <typename in_value_t, typename integral_t, typename out_value_t>
//     requires std::is_same_v<std::make_unsigned_t<std::decay_t<in_value_t>>, unsigned char> &&
//              std::is_same_v<std::make_unsigned_t<std::decay_t<out_value_t>>, unsigned char>
// inline void
// write(in_value_t * & iptr, integral_t const n, out_value_t * optr)
// {
//     std::memcpy(optr, iptr, n);
//     iptr += n;
// }

// write for more complex values (defer to write of iterator value)
// used for Strings of Pairs
//template <typename input_t, typename output_t, typename integral_t>
////inline SEQAN_FUNC_ENABLE_IF(Or<Is<OutputStreamConcept<TTarget> >, Is<ContainerConcept<TTarget> > >, void)
////inline SEQAN_FUNC_ENABLE_IF(And<
////                            Is<IntegerConcept<TSize> >,
////                            Not< Is<Convertible<typename Value<TTarget>::Type,
////                            typename Value<TFwdIterator>::Type> > > >, void)
//
//// not input::value_type not convertible with output_t::value_type
//inline void
//write(input_t  & i_iter,
//      output_t && o_iter,
//      integral_t n)
//{
//    for (; n > static_cast<TSize>(0); --n, ++i_iter)
//    {
//        put(*i_iter, std::forward<output_t>(o_iter));
//        put(' ', std::forward<output_t>(o_iter));
//    }
//}

template <typename input_t, typename integral_t, typename target_t>
inline void
write(input_t && input,
      integral_t const n,
      target_t & output)
    requires input_iterator_concept<std::remove_reference_t<input_t>> &&
             !iterator_concept<target_t> &&
             assignable_concept<value_type_t<target_t> &, value_type_t<std::remove_reference_t<input_t>>>
{
    write(std::forward<input_t>(input), n, output_iterator(output));
}

// ----------------------------------------------------------------------------
// Function write(TContainer) but not container of container
// ----------------------------------------------------------------------------

template <typename container_t, typename target_t>
inline void
write(container_t const & input,
      target_t & output)
    requires sized_range_concept<container_t> &&
             assignable_concept<value_type_t<target_t>&, value_type_t<container_t>>
{
    using std::size;
    write(std::get<0>(input_iterator(input)), size(input), output);
}

// ----------------------------------------------------------------------------
// Function read()
// ----------------------------------------------------------------------------

template <typename input_t,
          typename output_t,
          typename predicate_t>
//    requires output_iterator_concept<std::decay_t<target_t>> &&
//             input_range_concept<std::decay_t<input_t>> &&
//             predicate_concept<std::decay_t<predicate_t>>
inline void
get(input_t & curr_it,
    input_t const & end_it,
    output_t && output_it,
    predicate_t & check_func)
{
    if (/*SEQAN_UNLIKELY*/(curr_it == end_it))
        throw unexpected_end_error();

    AssertFunctor<predicate_t, parse_error> asserter(check_func);

    asserter(*curr_it);
    *output_it = *curr_it;
    ++curr_it;
}

template <typename input_t, typename output_t>
//    requires output_iterator_concept<std::decay_t<target_t>> &&
//             input_range_concept<std::decay_t<input_t>>
inline void
get(input_t & curr_it,
    input_t const & end_it,
    output_t && output_it)
{
    get(curr_it, end_it, std::forward<output_t>(output_it), always_true);
}

// ----------------------------------------------------------------------------
// Function _readUntil(); Element-wise
// ----------------------------------------------------------------------------

template <typename input_t,
          typename output_t,
          typename stop_predicate_t,
          typename ignore_predicate_t>
inline void
_get_until(input_t & i_iter,
           input_t const & i_end,
           output_t && output_it,
           stop_predicate_t && stop_func,
           ignore_predicate_t && ignore_func)
{
    typename std::remove_const<typename input_t::value_type>::type val;
    for (; i_iter != i_end; ++i_iter)
    {
        if (/*SEQAN_UNLIKELY(*/stop_predicate_t(val = *i_iter))/*)*/
            return;
        if (/*SEQAN_LIKELY(*/!ignoreFunctor(val))/*)*/
            put(val, std::forward<output_it>(output_it));
    }
}

// ----------------------------------------------------------------------------
// Function _readUntil(); Chunked
// ----------------------------------------------------------------------------
/*
 template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor, typename TIValue, typename TOValue>
 inline void _readUntil(TTarget &target,
 TFwdIterator &iter,
 TStopFunctor &stopFunctor,
 TIgnoreFunctor &ignoreFunctor,
 Range<TIValue*> *,
 Range<TOValue*> *)
 requires output_iterator_concept<TTarget> &&
 input_iterator_concept<TFwdIterator>
 {
 //
 Range<TOValue*> ochunk(NULL, NULL);
 TOValue* SEQAN_RESTRICT optr = NULL;

 Range<TIValue*> ichunk;
 for (; !atEnd(iter); )
 {
 getChunk(ichunk, iter, Input());  // get from gptr to egptr
 const TIValue* SEQAN_RESTRICT iptr = ichunk.begin;
 SEQAN_ASSERT(iptr < ichunk.end);

 for (; iptr != ichunk.end; ++iptr)
 {
 if (SEQAN_UNLIKELY(stopFunctor(*iptr)))
 {
 iter += iptr - ichunk.begin;               // advance input iterator
 advanceChunk(target, optr - ochunk.begin); // extend target string size
 return;
 }

 if (SEQAN_UNLIKELY(ignoreFunctor(*iptr)))
 continue;

 // construct values in reserved memory
 if (SEQAN_UNLIKELY(optr == ochunk.end))  // in case chunk is full.
 {
 // do a pbump or for container do a reserve.
 advanceChunk(target, optr - ochunk.begin);  // pbump or set
 // reserve memory for the worst-case
 // TODO(weese):Document worst-case behavior
 reserveChunk(target, length(ichunk), Output());  // reserves, but does not create.
 getChunk(ochunk, target, Output());
 optr = ochunk.begin;
 SEQAN_ASSERT(optr < ochunk.end);
 }
 //
 from_char(*optr++, *iptr);
 }
 iter += iptr - ichunk.begin;                       // advance input iterator
 }
 advanceChunk(target, optr - ochunk.begin);
 }
*/

// ----------------------------------------------------------------------------
// Function readUntil()
// ----------------------------------------------------------------------------

template <typename input_t,
          typename target_t,
          typename stop_predicate_t,
          typename ignore_predicate_t>
inline void
get_until(input_t & curr,
          input_t const & end,
          target_t && output_it,
          stop_predicate_t && stop_func,
          ignore_predicate_t && ignore_func)
{
    _read_until(curr, end, std::forward<target_t>(output_it),
                std::forward<stop_predicate_t>(stop_func), std::forward<ignore_predicate_t>(ignore_func));
}

// ----------------------------------------------------------------------------
// Function readUntil(); Not ignoring
// ----------------------------------------------------------------------------

template <typename input_t,
          typename target_t,
          typename stop_predicate_t>
inline void
get_until(input_t & curr,
          input_t const & end,
          target_t && output_it,
          stop_predicate_t && stop_func)
{
    read_until(curr, end, std::forward<target_t>(output_it),
               std::forward<stop_predicate_t>(stop_func), always_false);
}

// ----------------------------------------------------------------------------
// Function read_line()
// ----------------------------------------------------------------------------

template <typename input_t,
          typename target_t>
inline void
get_line(input_t & i_iter,
         input_t const & i_end,
         target_t && output_it)
{
    get_until(i_iter, i_end, std::forward<target_t>(output_it), IsNewline{});

    // consume "\r\n.", "\r[!\n]" or "\n."

    if (/*SEQAN_UNLIKELY(*/i_iter == i_end)/*)*/
        return;

    // If the current character is Line Feed ('\r') then this can be an ANSI or a Mac line ending.
    if (*i_iter == '\r')
    {
        ++i_iter;     // consume the found newline
        if (/*SEQAN_UNLIKELY(*/i_iter == i_end)/*)*/
            return;
    }

    // Unix Carriage Return ('\n') is the simplest case.
    if (*i_iter == '\n')
        ++i_iter;     // consume the found newline
}

// ----------------------------------------------------------------------------
// Function read()
// ----------------------------------------------------------------------------

template <typename input_t,
          typename output_t,
          typename integral_t,
          typename ignore_predicate_t>
inline void
read(input_t & i_iter,
     input_t const & i_end,
     output_t && o_iter,
     integral_t const n,
     ignore_predicate_t && ignore_func)
{
    // define counting stop_functor, which every time it's operator is called decrements n and checks if it is 0.
    auto stopper = [n](auto const &) { return n-- == 0; };
    get_until(i_iter, i_end, std::forward<output_t>(o_iter), stopper, std::forward<ignore_predicate_t>(ignore_func));
}

template <typename input_t,
          typename output_t,
          typename integral_t>
inline void
read(input_t & i_iter,
     input_t const & i_end,
     output_t && o_iter,
     integral_t const n)
{
    // define counting stop_functor, which every time it's operator is called decrements n and checks if it is 0.
    get(i_iter, i_end, std::forward<output_t>(o_iter), n, always_true);
}

// ----------------------------------------------------------------------------
// Function readRawByte()
// ----------------------------------------------------------------------------

template <typename input_t, typename value_t>
inline void
read_raw_pod(input_t & i_iter, value_t && value)
{
    // TODO(rrahn): need a back_insert_iterator adaptor for c-style array.
    write(i_iter, static_cast<char*>(&value), sizeof(value_t));
}

// ----------------------------------------------------------------------------
// Function _skipUntil(); Element-wise
// ----------------------------------------------------------------------------

template <typename input_t,
          typename stop_predicate_t>
inline void
_skip_until(input_t & i_iter,
            input_t const & i_end,
            stop_predicate_t && stop_func)
{
    for (; (i_iter != i_end) && !stop_func(*i_iter); ++i_iter)
    {}
}

// ----------------------------------------------------------------------------
// Function _skipUntil(); Chunked
// ----------------------------------------------------------------------------

/*
// TODO(rrahn): Replace Range<TValue*> * version
template <typename TFwdIterator, typename TStopFunctor, typename TValue>
inline void _skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor, Range<TValue*> *)
{
    // adapt the fwd_iterator interface.
    typedef typename Value<TFwdIterator>::Type TIValue;

    for (; !atEnd(iter); )
    {
        // We operate on bare pointers.
        // can we also have a view?
        Range<TIValue const *> ichunk;
        // TODO(rrahn): getChunk -> implements the chunking interface, what should returned here?

        // this returns the chunk based on the iterator.
        // This means the iter needs access to it's underlying container.

        // We get the remaining characters from the buffer
        // We can model this with diff betwenn pptr() and egptr()
        getChunk(ichunk, iter, Input());

        // The rest of the buffer.
        // for vector<char> -> this is the remaining part from it to end of iterator.
        // we need to model this through the iterator.

        // Chunk on iterator gives back a
        SEQAN_ASSERT(!empty(ichunk));

        const TIValue* SEQAN_RESTRICT ptr = ichunk.begin;

        for (; ptr != ichunk.end; ++ptr)
        {
            // Now we can run the pointer here.
            if (SEQAN_UNLIKELY(stopFunctor(*ptr)))
            {
                iter += ptr - ichunk.begin;    // advance input iterator
                return;
            }
        }

        iter += ptr - ichunk.begin;            // advance input iterator
    }
}
*/
// ----------------------------------------------------------------------------
// Function skipUntil()
// ----------------------------------------------------------------------------

template <typename input_t,
          typename stop_predicate_t>
inline void
skip_until(input_t & i_iter,
           input_t const & i_end,
           stop_predicate_t && stop_func)
{
    _skip_until(i_iter, i_end, std::forward<stop_predicate_t>(stop_func));
}

// ----------------------------------------------------------------------------
// Function skip
// ----------------------------------------------------------------------------

template <typename input_t,
          typename unexpected_predicate_t>
inline void
skip(input_t & i_iter,
     input_t const & i_end,
     unexpected_predicate_t && unexpected_func)
{
    AssertFunctor<unexpected_predicate_t, parse_error> asserter{unexpected_func};

    if (/*SEQAN_UNLIKELY*/(i_iter == i_end))
        throw unexpected_end_error{};

    asserter(*i_iter);
    ++i_iter;
}

template <typename input_t>
inline void
skip(input_t & i_iter,
     input_t const & i_end)
{
    skip(i_iter, i_end, always_true);
}

// ----------------------------------------------------------------------------
// Function skipLine()
// ----------------------------------------------------------------------------

template <typename input_t>
inline void
skip_line(input_t & i_iter,
          input_t const & i_end)
{
    skip_until(i_iter, i_end, IsNewline());

    // consume "\r\n.", "\r[!\n]" or "\n."

    if (/*SEQAN_UNLIKELY(*/i_iter == i_end)/*)*/
        return;

    // If the current character is Line Feed ('\r') then this can be an ANSI or a Mac line ending.
    if (*i_iter == '\r')
    {
        ++i_iter;     // consume the found newline
        if (/*SEQAN_UNLIKELY(*/i_iter == i_end)/*)*/
            return;
    }

    // Unix Carriage Return ('\n') is the simplest case.
    if (*i_iter == '\n')
        ++i_iter;     // consume the found newline
}

// ----------------------------------------------------------------------------
// Function writeWrappedString()
// ----------------------------------------------------------------------------

//template <typename input_range_t, typename output_t, typename integral_t>
//inline void
//write_wrapped(input_range_t & in,
//              output_t && o_iter,
//              integral_t const line_length)
//{
//    // TODO(rrahn): Global Metafunction.
//    using size_type = input_range_t::size_type;
//
//    auto i_iter = in_begin(in);
//    size_type chars_left = std::size(in);
//    size_type chars_per_line;
//    size_type line_length_ = (line_length == 0)? maxValue<size_type>() : line_length;
//    do
//    {
//        chars_per_line = std::min(chars_left, line_length_);
//        write(i_iter, std::forward<output_t>(o_iter), chars_per_line);
//        put('\n', std::forward<output_t>(o_iter));
//        chars_left -= chars_per_line;
//    }
//    while (chars_left != 0);
//}

//// ----------------------------------------------------------------------------
//// Function findFirst()
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TFunctor>
//inline typename Position<TContainer>::Type
//findFirst(TContainer const &cont, TFunctor const &func)
//{
//    typename Iterator<TContainer const, Rooted>::Type iter = begin(cont, Rooted());
//    skipUntil(iter, func);
//    return iter - begin(cont, Rooted());
//}
//
//template <typename TContainer>
//inline typename Position<TContainer>::Type
//findFirst(TContainer const &cont, typename Value<TContainer>::Type const &val)
//{
//    EqualsDynamicValue<typename Value<TContainer>::Type> func(val);
//    return findFirst(cont, func);
//}
//
//// ----------------------------------------------------------------------------
//// Function findLast()
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TFunctor>
//inline typename Position<TContainer>::Type
//findLast(TContainer const &cont, TFunctor const &func)
//{
//    typedef ModifiedString<TContainer const, ModReverse> TRevContainer;
//
//    SEQAN_CONCEPT_ASSERT((IntegerConcept<typename Position<TContainer>::Type>));
//
//    // search from back to front
//    TRevContainer rev(cont);
//    typename Iterator<TRevContainer, Rooted>::Type iter = begin(rev, Rooted());
//    skipUntil(iter, func);
//
//    if (atEnd(iter))
//        return -1;
//
//    return host(iter) - begin(cont, Rooted());
//}
//
//template <typename TContainer>
//inline typename Position<TContainer>::Type
//findLast(TContainer const &cont, typename Value<TContainer>::Type const &val)
//{
//    EqualsDynamicValue<typename Value<TContainer>::Type> func(val);
//    return findLast(cont, func);
//}
//
//// ----------------------------------------------------------------------------
//// Function cropAfterFirst(); crop after first occurrence (including it)
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TFunctor>
//inline void
//cropAfterFirst(TContainer &cont, TFunctor const &func)
//{
//    resize(cont, findFirst(cont, func));
//}
//
//// ----------------------------------------------------------------------------
//// Function cropAfterLast(); crop after last occurrence (excluding it)
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TFunctor>
//inline void
//cropAfterLast(TContainer &cont, TFunctor const &func)
//{
//    resize(cont, findLast(cont, func) + 1);
//}
//
//// ----------------------------------------------------------------------------
//// Function cropBeforeFirst(); crop before first occurrence (excluding it)
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TFunctor>
//inline void
//cropBeforeFirst(TContainer &cont, TFunctor const &func)
//{
//    erase(cont, 0, findFirst(cont, func));
//}
//
//// ----------------------------------------------------------------------------
//// Function cropBeforeLast(); crop before first occurrence (including it)
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TFunctor>
//inline void
//cropBeforeLast(TContainer &cont, TFunctor const &func)
//{
//    erase(cont, 0, findLast(cont, func) + 1);
//}
//// ----------------------------------------------------------------------------
//// Function cropOuter(); crop after last occurrence (excluding it)
//// ----------------------------------------------------------------------------
//
//template <typename TContainer, typename TFunctor>
//inline void
//cropOuter(TContainer &cont, TFunctor const &func)
//{
//    cropAfterLast(cont, NotFunctor<TFunctor>(func));
//    cropBeforeFirst(cont, NotFunctor<TFunctor>(func));
//}

// --------------------------------------------------------------------------
// Function strSplit()
// --------------------------------------------------------------------------

/*!
 * @fn StringSet#strSplit
 * @brief Split a sequence using a delimiter and append the results to a target string set
 *
 * @signature void strSplit(result, sequence[, sep[, allowEmptyStrings[, maxSplit]]]);
 *
 * @param[out] result           The resulting string set (can be any ContainerOfContainer, also STL)
 * @param[in]  sequence         The sequence to split.
 * @param[in]  sep              The splitter to use (default <tt>' '</tt>).
 * @param[in]  allowEmptyString Whether or not to allow empty strings (<tt>bcd gitool</tt>, defaults to <tt>true</tt> iff
 *                              <tt>sep</tt> is given).
 * @param[in]  maxSplit         The maximal number of split operations to do if given.
 */

//template <typename sequence_t, typename TFunctor, typename TSize>
//    requires forward_range_concept<sequence_t> &&
//
//inline auto
//split_by(TResult & result,
//         TSequence const & sequence,
//         TFunctor const & sep,
//         bool const allowEmptyStrings,
//         TSize maxSplit)
//{
//    typedef typename Iterator<TSequence const, Standard>::Type TIter;
//    typedef std::conditional_t<Is<StlContainerConcept<TResult>>::VALUE,
//    TSequence,
//    decltype(infix(sequence, 0, 1))> TResultValue;
//
//    TIter itBeg = begin(sequence, Standard());
//    TIter itEnd = end(sequence, Standard());
//    TIter itFrom = itBeg;
//
//    if (maxSplit == 0)
//    {
//        appendValue(result, sequence);
//        return;
//    }
//
//    for (TIter it = itBeg; it != itEnd; ++it)
//        if (sep(getValue(it)))
//        {
//            if (allowEmptyStrings || itFrom != it)
//            {
//                appendValue(result, static_cast<TResultValue>(infix(sequence, itFrom - itBeg, it - itBeg)));
//                if (--maxSplit == 0)
//                {
//                    if (!allowEmptyStrings)
//                    {
//                        while (it != itEnd && sep(getValue(it)))
//                            ++it;
//                    }
//                    else
//                        ++it;
//
//                    if (it != itEnd)
//                        appendValue(result, static_cast<TResultValue>(infix(sequence, itFrom - itBeg, it - itBeg)));
//
//                        return;
//                }
//            }
//            itFrom = it + 1;
//        }
//
//    if (allowEmptyStrings || itFrom != itEnd)
//        appendValue(result, static_cast<TResultValue>(infix(sequence, itFrom - itBeg, itEnd - itBeg)));
//}
//
//template <typename TResult, typename TSequence, typename TFunctor>
//inline SEQAN_FUNC_ENABLE_IF(And<Is<ContainerConcept<TResult> >,
//                            Is<ContainerConcept<typename Value<TResult>::Type > > >, void)
//strSplit(TResult & result, TSequence const & sequence, TFunctor const & sep, bool const allowEmptyStrings)
//{
//    strSplit(result, sequence, sep, allowEmptyStrings, maxValue<typename Size<TSequence>::Type>());
//}
//
//template <typename TResult, typename TSequence, typename TFunctor>
//inline SEQAN_FUNC_ENABLE_IF(And<Is<ContainerConcept<TResult> >,
//                            Is<ContainerConcept<typename Value<TResult>::Type > > >, void)
//strSplit(TResult & result, TSequence const & sequence, TFunctor const & sep)
//{
//    strSplit(result, sequence, sep, true);
//}
//
//template <typename TResult, typename TSequence>
//inline SEQAN_FUNC_ENABLE_IF(And<Is<ContainerConcept<TResult> >,
//                            Is<ContainerConcept<typename Value<TResult>::Type > > >, void)
//strSplit(TResult & result, TSequence const & sequence)
//{
//    strSplit(result, sequence, EqualsChar<' '>(), false);
//}
}  // namespace seqan3::detail
