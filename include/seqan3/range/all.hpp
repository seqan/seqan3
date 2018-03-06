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

/*!\file
 * \brief Meta-header for the \link range range module \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/concept.hpp>
#include <seqan3/range/action/all.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/range/view/all.hpp>

/*!\defgroup range Range
 * \brief The range module contains containers, views and actions.
 *
 * *Ranges* are an abstraction of "a collection of items", or "something iterable". The most basic definition
 * requires only the existence of begin() and end() on the range. See range/concept.hpp for the different range
 * concepts.
 *
 * *Containers* are ranges that own their elements. SeqAn3 makes use of standard STL containers like std::vector,
 * but also implements some custom containers. See range/container.hpp for more details.
 *
 * *Views* are "lazy range combinators" that provide operations on other ranges, e.g. containers, but do so on-demand,
 * i.e. views don't own elements, but return (mutated) elements on request. This is similar to how iterators can
 * provide different behaviours on the same underlying data structure (while not actually changing it). See
 * range/view.hpp for more details.
 *
 * *Actions* on the other hand are "eager range combinators", i.e. they immediately change the underlying range
 * they are applied to. See range/action.hpp for more details.
 *
 * Both views and actions can be chained via the pipe operator.
 *
 * \sa range.hpp
 * \sa https://ericniebler.github.io/range-v3/index.html
 */
