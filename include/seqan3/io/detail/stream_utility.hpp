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
#include <string>
#include <tuple>
#include <stdexcept>

#include "seqan3/alphabet/alphabet.hpp"

namespace seqan3::detail {

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


using  RuntimeError = std::runtime_error;

struct ParseError : RuntimeError
{
    template <typename TString>
    ParseError(TString const & message) :
    RuntimeError(message)
    {}
};

struct UnexpectedEnd : ParseError
{

    UnexpectedEnd() :
    ParseError("Unexpected end of input.")
    {}
};

// ----------------------------------------------------------------------------
// Exception EmptyFieldError
// ----------------------------------------------------------------------------

struct EmptyFieldError : ParseError
{
    EmptyFieldError(std::string fieldName):
    ParseError(fieldName + " field was empty.")
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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _skipUntil(); Element-wise
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TStopFunctor, typename TChunk>
inline void _skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor, TChunk const &)
{
    for (; !atEnd(iter) && !stopFunctor(*iter); ++iter) ;
}

// ----------------------------------------------------------------------------
// Function _skipUntil(); Chunked
// ----------------------------------------------------------------------------

/*
 template <typename TContainer, typename TValue, typename TSpec>
 struct Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >
 {
 typedef Range<TValue*> Type;
 };
 
 template <typename TIterator, typename TSize>
 inline void advanceChunk(TIterator &iter, TSize size)
 {
 iter += size;
 }
 
 
 template <typename TContainer, typename TSpec, typename TSize>
 inline void reserveChunk(Iter<TContainer, TSpec> &iter, TSize size, Output)
 {
 typedef Iter<TContainer, TSpec> TIter;

 TContainer &cont = container(iter);
 typename Position<TIter>::Type pos = position(iter);
 typename Size<TIter>::Type len = length(cont);

 if (pos < len)
 return;

 len += size;
 if (len <= capacity(cont))
 return;

 reserve(cont, len);
 setPosition(iter, pos);
 }

 template <typename TContainer, typename TSpec, typename TSize>
 inline void reserveChunk(Iter<TContainer, TSpec> &, TSize, Input)
 {}


 template <typename TChunk, typename TContainer, typename TValue, typename TSpec>
 inline void
 getChunk(TChunk &result, Iter<TContainer, AdaptorIterator<TValue*, TSpec> > &rootedIter, Input)
 {
 return assignRange(result, hostIterator(rootedIter), end(container(rootedIter), Standard()));
 }

 template <typename TChunk, typename TContainer, typename TValue, typename TSpec>
 inline void
 getChunk(TChunk &result, Iter<TContainer, AdaptorIterator<TValue*, TSpec> > &rootedIter, Output)
 {
 TContainer &cont = container(rootedIter);
 return assignRange(result, hostIterator(rootedIter), begin(cont, Standard()) + capacity(cont));
 }
*/

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

// ----------------------------------------------------------------------------
// Function skipUntil()
// ----------------------------------------------------------------------------

// I need to wrap the
template <typename TFwdIterator, typename TStopFunctor>
inline void skipUntil(TFwdIterator &iter, TStopFunctor &stopFunctor)
{
//    typedef typename Chunk<TFwdIterator>::Type* TIChunk;

    _skipUntil(iter, stopFunctor, std::ignore);
}

template <typename TFwdIterator, typename TStopFunctor>
inline void skipUntil(TFwdIterator &iter, TStopFunctor const &stopFunctor)
{
//    typedef typename Chunk<TFwdIterator>::Type* TIChunk;

    TStopFunctor stopFunctor_ = stopFunctor;
    _skipUntil(iter, stopFunctor_, std::ignore);
}

// ----------------------------------------------------------------------------
// Function skipOne()
// ----------------------------------------------------------------------------

template <typename TFwdIterator, typename TFunctor>
inline void skipOne(TFwdIterator &iter, TFunctor &functor)
{
    AssertFunctor<TFunctor, ParseError> asserter(functor);

    if (SEQAN_UNLIKELY(atEnd(iter)))
        throw UnexpectedEnd();

    asserter(*iter);
    ++iter;
}

template <typename TFwdIterator, typename TFunctor>
inline void skipOne(TFwdIterator &iter, TFunctor const &functor)
{
    TFunctor func(functor);
    skipOne(iter, func);
}

template <typename TFwdIterator>
inline void skipOne(TFwdIterator &iter)
{
    skipOne(iter, [](auto const &){ return true; });
}

// ----------------------------------------------------------------------------
// Function _readUntil(); Element-wise
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor, typename TIChunk, typename TOChunk>
inline void
_readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor, TIgnoreFunctor &ignoreFunctor, TIChunk, TOChunk)
{
    typename std::remove_const<typename Value<TFwdIterator>::Type>::type val;
    for (; !atEnd(iter); ++iter)
    {
        if (/*SEQAN_UNLIKELY(*/stopFunctor(val = *iter))/*)*/
            return;
        if (/*SEQAN_LIKELY(*/!ignoreFunctor(val))/*)*/
            writeValue(target, val);
    }
}

// ----------------------------------------------------------------------------
// Function _readUntil(); Chunked
// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------
// Function readUntil()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor>
inline void readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor, TIgnoreFunctor &ignoreFunctor)
{
//    typedef typename Chunk<TFwdIterator>::Type*         TIChunk;
//    typedef typename Iterator<TTarget, Rooted>::Type    TTargetIter;
//    typedef typename Chunk<TTargetIter>::Type*          TOChunk;

    _readUntil(target, iter, stopFunctor, ignoreFunctor, std::ignore, std::ignore);
}

template <typename TTarget, typename TFwdIterator, typename TStopFunctor, typename TIgnoreFunctor>
inline void
readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor const &stopFunctor, TIgnoreFunctor const &ignoreFunctor)
{
    TStopFunctor stopFunctor_ = stopFunctor;
    TIgnoreFunctor ignoreFunctor_ = ignoreFunctor;
    readUntil(target, iter, stopFunctor_, ignoreFunctor_);
}

// ----------------------------------------------------------------------------
// Function readUntil(); Not ignoring
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TStopFunctor>
inline void readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor const &stopFunctor)
{
    readUntil(target, iter, stopFunctor, False());
}

template <typename TTarget, typename TFwdIterator, typename TStopFunctor>
inline void readUntil(TTarget &target, TFwdIterator &iter, TStopFunctor &stopFunctor)
{
    readUntil(target, iter, stopFunctor, [](auto const &){ return false; });
}

// ----------------------------------------------------------------------------
// Function readOne()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator, typename TFunctor>
inline void readOne(TTarget & target, TFwdIterator &iter, TFunctor &functor)
{
    if (/*SEQAN_UNLIKELY(*/atEnd(iter))/*)*/
        throw UnexpectedEnd();

    AssertFunctor<TFunctor, ParseError> asserter(functor);

    asserter(*iter);
    target = *iter;
    ++iter;
}

template <typename TTarget, typename TFwdIterator, typename TFunctor>
inline void readOne(TTarget & target, TFwdIterator &iter, TFunctor const &functor)
{
    TFunctor func(functor);
    readOne(target, iter, func);
}

template <typename TTarget, typename TFwdIterator>
inline void readOne(TTarget & target, TFwdIterator &iter)
{
    if (/*SEQAN_UNLIKELY(*/atEnd(iter))/*)*/
        throw UnexpectedEnd();

    target = *iter;
    ++iter;
}

// ----------------------------------------------------------------------------
// Function readRawByte()
// ----------------------------------------------------------------------------

//TODO(singer) to be revised
template <typename TValue, typename TFwdIterator>
inline void readRawPod(TValue & value, TFwdIterator &srcIter)
{
    write((char*)&value, srcIter, sizeof(TValue));
}

// ----------------------------------------------------------------------------
// Function readLine()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TFwdIterator>
inline void readLine(TTarget &target, TFwdIterator &iter)
{
    readUntil(target, iter, IsNewline());

    // consume "\r\n.", "\r[!\n]" or "\n."

    if (/*SEQAN_UNLIKELY(*/atEnd(iter))/*)*/
        return;

    // If the current character is Line Feed ('\r') then this can be an ANSI or a Mac line ending.
    if (*iter == '\r')
    {
        ++iter;     // consume the found newline
        if (/*SEQAN_UNLIKELY(*/atEnd(iter))/*)*/
            return;
    }

    // Unix Carriage Return ('\n') is the simplest case.
    if (*iter == '\n')
        ++iter;     // consume the found newline
}

// ----------------------------------------------------------------------------
// Function skipLine()
// ----------------------------------------------------------------------------

template <typename TFwdIterator>
inline void skipLine(TFwdIterator &iter)
{
    skipUntil(iter, IsNewline());

    // consume "\r\n.", "\r[!\n]" or "\n."

    if (/*SEQAN_UNLIKELY(*/atEnd(iter))/*)*/
        return;

    // If the current character is Line Feed ('\r') then this can be an ANSI or a Mac line ending.
    if (*iter == '\r')
    {
        ++iter;     // consume the found newline
        if (/*SEQAN_UNLIKELY(*/atEnd(iter))/*)*/
            return;
    }

    // Unix Carriage Return ('\n') is the simplest case.
    if (*iter == '\n')
        ++iter;     // consume the found newline
}

// ----------------------------------------------------------------------------
// Function writeWrappedString()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSequence, typename TSize>
inline void writeWrappedString(TTarget & target, TSequence const & seq, TSize lineLength)
{
    typedef typename TSequence::size_type TSeqSize;
    typedef typename Iterator<TSequence const, Rooted>::Type TIter;

    TIter iter = begin(seq, Rooted());
    TSeqSize charsLeft = length(seq);
    TSeqSize charsPerLine;
    TSeqSize lineLength_ = (lineLength == 0)? maxValue<TSeqSize>() : lineLength;

    do
    {
        charsPerLine = std::min(charsLeft, lineLength_);
        write(target, iter, charsPerLine);
        writeValue(target, '\n');
        charsLeft -= charsPerLine;
    }
    while (charsLeft != 0);
}

// ----------------------------------------------------------------------------
// Function findFirst()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline typename Position<TContainer>::Type
findFirst(TContainer const &cont, TFunctor const &func)
{
    typename Iterator<TContainer const, Rooted>::Type iter = begin(cont, Rooted());
    skipUntil(iter, func);
    return iter - begin(cont, Rooted());
}

template <typename TContainer>
inline typename Position<TContainer>::Type
findFirst(TContainer const &cont, typename Value<TContainer>::Type const &val)
{
    EqualsDynamicValue<typename Value<TContainer>::Type> func(val);
    return findFirst(cont, func);
}

// ----------------------------------------------------------------------------
// Function findLast()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline typename Position<TContainer>::Type
findLast(TContainer const &cont, TFunctor const &func)
{
    typedef ModifiedString<TContainer const, ModReverse> TRevContainer;

    SEQAN_CONCEPT_ASSERT((IntegerConcept<typename Position<TContainer>::Type>));

    // search from back to front
    TRevContainer rev(cont);
    typename Iterator<TRevContainer, Rooted>::Type iter = begin(rev, Rooted());
    skipUntil(iter, func);

    if (atEnd(iter))
        return -1;

    return host(iter) - begin(cont, Rooted());
}

template <typename TContainer>
inline typename Position<TContainer>::Type
findLast(TContainer const &cont, typename Value<TContainer>::Type const &val)
{
    EqualsDynamicValue<typename Value<TContainer>::Type> func(val);
    return findLast(cont, func);
}

// ----------------------------------------------------------------------------
// Function cropAfterFirst(); crop after first occurrence (including it)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline void
cropAfterFirst(TContainer &cont, TFunctor const &func)
{
    resize(cont, findFirst(cont, func));
}

// ----------------------------------------------------------------------------
// Function cropAfterLast(); crop after last occurrence (excluding it)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline void
cropAfterLast(TContainer &cont, TFunctor const &func)
{
    resize(cont, findLast(cont, func) + 1);
}

// ----------------------------------------------------------------------------
// Function cropBeforeFirst(); crop before first occurrence (excluding it)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline void
cropBeforeFirst(TContainer &cont, TFunctor const &func)
{
    erase(cont, 0, findFirst(cont, func));
}

// ----------------------------------------------------------------------------
// Function cropBeforeLast(); crop before first occurrence (including it)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline void
cropBeforeLast(TContainer &cont, TFunctor const &func)
{
    erase(cont, 0, findLast(cont, func) + 1);
}
// ----------------------------------------------------------------------------
// Function cropOuter(); crop after last occurrence (excluding it)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline void
cropOuter(TContainer &cont, TFunctor const &func)
{
    cropAfterLast(cont, NotFunctor<TFunctor>(func));
    cropBeforeFirst(cont, NotFunctor<TFunctor>(func));
}

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
 * @param[in]  allowEmptyString Whether or not to allow empty strings (<tt>bool</tt>, defaults to <tt>true</tt> iff
 *                              <tt>sep</tt> is given).
 * @param[in]  maxSplit         The maximal number of split operations to do if given.
 */

template <typename TResult, typename TSequence, typename TFunctor, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(And<Is<ContainerConcept<TResult> >,
                            Is<ContainerConcept<typename Value<TResult>::Type > > >, void)
strSplit(TResult & result,
         TSequence const & sequence,
         TFunctor const & sep,
         bool const allowEmptyStrings,
         TSize maxSplit)
{
    typedef typename Iterator<TSequence const, Standard>::Type TIter;
    typedef std::conditional_t<Is<StlContainerConcept<TResult>>::VALUE,
    TSequence,
    decltype(infix(sequence, 0, 1))> TResultValue;
    
    TIter itBeg = begin(sequence, Standard());
    TIter itEnd = end(sequence, Standard());
    TIter itFrom = itBeg;
    
    if (maxSplit == 0)
    {
        appendValue(result, sequence);
        return;
    }
    
    for (TIter it = itBeg; it != itEnd; ++it)
        if (sep(getValue(it)))
        {
            if (allowEmptyStrings || itFrom != it)
            {
                appendValue(result, static_cast<TResultValue>(infix(sequence, itFrom - itBeg, it - itBeg)));
                if (--maxSplit == 0)
                {
                    if (!allowEmptyStrings)
                    {
                        while (it != itEnd && sep(getValue(it)))
                            ++it;
                    }
                    else
                        ++it;
                    
                    if (it != itEnd)
                        appendValue(result, static_cast<TResultValue>(infix(sequence, itFrom - itBeg, it - itBeg)));
                        
                        return;
                }
            }
            itFrom = it + 1;
        }
    
    if (allowEmptyStrings || itFrom != itEnd)
        appendValue(result, static_cast<TResultValue>(infix(sequence, itFrom - itBeg, itEnd - itBeg)));
}

template <typename TResult, typename TSequence, typename TFunctor>
inline SEQAN_FUNC_ENABLE_IF(And<Is<ContainerConcept<TResult> >,
                            Is<ContainerConcept<typename Value<TResult>::Type > > >, void)
strSplit(TResult & result, TSequence const & sequence, TFunctor const & sep, bool const allowEmptyStrings)
{
    strSplit(result, sequence, sep, allowEmptyStrings, maxValue<typename Size<TSequence>::Type>());
}

template <typename TResult, typename TSequence, typename TFunctor>
inline SEQAN_FUNC_ENABLE_IF(And<Is<ContainerConcept<TResult> >,
                            Is<ContainerConcept<typename Value<TResult>::Type > > >, void)
strSplit(TResult & result, TSequence const & sequence, TFunctor const & sep)
{
    strSplit(result, sequence, sep, true);
}

template <typename TResult, typename TSequence>
inline SEQAN_FUNC_ENABLE_IF(And<Is<ContainerConcept<TResult> >,
                            Is<ContainerConcept<typename Value<TResult>::Type > > >, void)
strSplit(TResult & result, TSequence const & sequence)
{
    strSplit(result, sequence, EqualsChar<' '>(), false);
}
}  // namespace detail
}  // namespace seqan3
