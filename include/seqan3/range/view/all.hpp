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
 * \ingroup view
 * \brief Meta-header for the \link view view submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/concept.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/rank_to.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/range/view/trim.hpp>

/*!\defgroup view View
 * \brief Views are "lazy range combinators" that offer modified views onto other ranges.
 * \ingroup range
 *
 * \details
 *
 * SeqAn3 makes heavy use of views as defined in the
 * [Ranges Technical Specification](http://en.cppreference.com/w/cpp/experimental/ranges). From the original
 * documentation:  <i>"A view is a lightweight wrapper that presents a view of an underlying sequence of elements in
 * some custom way without mutating or copying it. Views are cheap to create and copy, and have non-owning reference
 * semantics. [...] The big advantage of ranges over iterators is their composability. They permit a functional style
 * of programming where data is manipulated by passing it through a series of combinators. In addition, the combinators
 * can be lazy, only doing work when the answer is requested, and purely functional, without mutating the original
 * data. This makes it easier to reason about your code, especially when writing concurrent programs."</i>
 *
 * The (less flexible) equivalent in SeqAn2 was the `ModifiedString<>`.
 *
 * Most views provided by SeqAn3 are specific to biological operations, like seqan3::view::trim which trims
 * sequences based on the quality or seqan3::view::complement which generates the complement of a nucleotide sequence.
 * But SeqAn3 also provides some general purpose views.
 *
 * \par Namespaces
 *
 *   * [All views from the range-v3 libary](https://ericniebler.github.io/range-v3/index.html#range-views) are available
 * in the namespace `ranges::view`.
 *
 *   * All SeqAn views are available in the namespace `seqan3::view`.
 *
 * \par Example
 *
 * Functional and pipe notations:
 * ```cpp
 * dna4_vector vec{"ACGGTC"_dna4};
 *
 * // these are synonymous:
 * auto vec_view1 = vec | view::complement;
 * auto vec_view2 = view::complement(vec);
 *
 * // both views "behave" like a collection of the elements 'T', 'G', 'C', 'C', 'A', 'G'
 * // but can be copied cheaply et cetera
 * ```
 *
 * Re-transform into a distinct container:
 * ```cpp
 * // just re-assign to a container
 * dna4_vector complemented = vec_view2;
 * assert(complemented == "TGCCAG"_dna4);
 *
 * // or immediately create on container
 * dna4_vector reversed = vec | ranges::view::reverse;
 * assert(complemented == "CTGGCA"_dna4);
 * ```
 *
 * Composability:
 * ```cpp
 * // views can be composed iteratively
 * auto vec_view3 = vec | ranges::view::reverse;
 * auto vec_view4 = vec_view3 | view::complement;
 *
 * // or in one line similar to the unix command line
 * auto vec_view5 = vec | view::complement | ranges::view::reverse;
 *
 * // vec_view4 and vec_view5 are the reverse complement of "ACGGTC": "GACCGT"
 * ```
 *
 * \par View properties
 *
 * All of these are documented for SeqAn3's views individually, but only explained here in detail:
 *
 * **Source-only views:** Most views operate on an input range and return a (modified) range, i.e. they can be placed
 * at the beginning, middle or end of a "pipe" of view operations. However, some views are limited to being at
 * the front ("source"), e.g. `ranges::view::single`, `ranges::view::concat` and `ranges::view::ints`.
 *
 * **Sink-only views:** The opposite of a *source-only view*. It can only be placed at the end of a pipe, i.e.
 * it operates on views, but does not actually return a view.
 *
 * **Deep views:** Some views are declared as "deeps views". This means, that in case they are given a range-of-range
 * as input (as opposed to just a range), they will apply their transformation on the innermost range (instead of
 * the outermost range which would be default). This is handy especially for alphabet-based transformations that you
 * wish to apply to a collection of sequences.
 *
 * **Input range requirements:** All views that are not *source-only views* make certain assumptions about their input.
 * The most basic assumption is that the range satisfies `seqan3::input_range_concept`, but some views require
 * stronger properties, e.g. `seqan3::random_access_range_concept`. *Note that these being* requirements *means that
 * they are the minimal set of properties assumed. Views may very well make use of stronger properties if available.*
 *
 * **Return range guarantees:** All views that are not *sink-only views* return a range that meets at least
 * `seqan3::input_range_concept` and of course also `seqan3::view_concept`. Most views also preserve stronger
 * properties, e.g. `seqan3::random_access_range_concept`, but this depends on the view. Some views also add
 * proporties not present on the input range, e.g. the range returned by `ranges::view::take_exactly` meets
 * `seqan3::sized_range_concept`, independent of whether this was met by the input range.
 *
 * **Input range's reference type:** The reference type is the type the elements of the input range are accessed by
 * (since dereferencing an iterator or calling operator[] returns the reference type). The reference type may or may
 * not actually contain a `&` (see below). For many SeqAn specific views additional concept requirements are defined
 * for the input range's reference type, e.g. seqan3::view::complement can only operate on ranges whose elements are
 * nucleotides (meet seqan3::nucleotide_concept). In some case the type may even be a specific type or the result
 * of a metafunction.
 *
 * **Return range's reference type:** Conversely certain views make guarantees on the concepts satisfied by the
 * return range's reference type or even always have a fixed type, e.g. seqan3::view::complement operates on
 * nucleotides and of course also returns nucleotides. However, and this is important to note, the reference type
 * of seqan3::view::complement has any actual `&` removed from the input ranges' reference type (if originally present).
 * This is because *new elements* are being generated. Other views like `ranges::view::reverse` also preserve the
 * `&` (if originally present), because the elements in the return view still point to the elements in the original
 * range (just in different order). This has the effect that through some combinations of views you can modify the
 * elements in the original range (if all views in the pipe preserve the reference type fully), but through others
 * you can't.
 *
 * \sa https://ericniebler.github.io/range-v3/index.html#range-views
 */

/*!
 * \namespace seqan3::view
 * \brief The SeqAn3 namespace for views.
 *
 * Since views often have name clashes with regular functions and ranges they are implemented in the sub
 * namespace `view`.
 *
 * See the \link view view submodule \endlink of the range module for more details.
 */
