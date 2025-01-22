// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link utility_views Utility / Views submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

/*!\defgroup utility_views Views
 * \brief Views are "lazy range combinators" that offer modified views onto other ranges.
 * \ingroup utility
 *
 * \details
 *
 * SeqAn makes heavy use of views as defined in the
 * [Ranges Technical Specification](https://en.cppreference.com/w/cpp/experimental/ranges). From the original
 * documentation:  <i>"A view is a lightweight wrapper that presents a view of an underlying sequence of elements in
 * some custom way without mutating or copying it. Views are cheap to create and copy, and have non-owning reference
 * semantics. [...] The big advantage of ranges over iterators is their composability. They permit a functional style
 * of programming where data is manipulated by passing it through a series of combinators. In addition, the combinators
 * can be lazy, only doing work when the answer is requested, and purely functional, without mutating the original
 * data. This makes it easier to reason about your code, especially when writing concurrent programs."</i>
 *
 * See the \link utility_range range module \endlink for how views relate to containers and decorators.
 *
 * Most views provided by SeqAn are specific to biological operations, like seqan3::views::trim_quality which trims
 * sequences based on the quality or seqan3::views::complement which generates the complement of a nucleotide sequence.
 * But SeqAn also provides some general purpose views.
 *
 * ### Namespaces
 *
 *   * [All views from the std libary](https://en.cppreference.com/w/cpp/ranges#Range_adaptors) are available
 * in the namespace `std::ranges::views` or in the namespace alias `std::views`.
 *
 *   * All SeqAn views are available in the namespace `seqan3::views`.
 *
 * ### Example
 *
 * Functional and pipe notations:
 * \include test/snippet/range/views/range_view_all_notation.cpp
 *
 * Re-transform into a distinct container:
 * \include test/snippet/range/views/range_view_all_retransform.cpp
 *
 * Composability:
 * \include test/snippet/range/views/range_view_all_composability.cpp
 *
 * ### Views vs view adaptors
 *
 * When talking about views, two different entities are often conflated:
 *
 *   1. the view (this is the type that is a range and meets std::ranges::view; it is what we refer to with
 * `auto vec_view` above)
 *   2. the view adaptor (this is the functor that returns the actual view based on it's parameters, including the
 * underlying range; in the above example `std::views::reverse` and `seqan3::views::complement` are view adaptors)
 *
 * The view adaptor also facilitates the piping behaviour. It is the only entity that is publicly documented and
 * the actual view type (the range type returned by the adaptor) is considered implementation defined.
 * The *properties* of the returned range are however specified and documented as part of the adaptor, see below.
 *
 * <sub>
 * An exception to this rule are views that don't work on an underlying range and can only be
 * placed at the beginning of a pipe of operations; they do not need an adaptor, because their constructor is
 * sufficient. This is not relevant for the documentation, though, we always document `views::foo`, independent of
 * whether `views::foo` is the adaptor that returns the "foo type" or whether `views::foo` is the "foo type".
 * </sub>
 *
 * ### View properties
 *
 * There are three view properties that are documented for a view, **only if that view fulfils them:**
 *
 * **Source-only views:** Most views operate on an underlying range and return a (modified) range, i.e. they can be placed
 * at the beginning, middle or end of a "pipe" of view operations. However, some views are limited to being at
 * the front ("source"), e.g. `std::views::empty`, `std::view::single` and `std::views::iota`. These views
 * are marked as "source-only" and have no `urng_t` column in the second table. The C++ standard calls these [Range
 * factories](https://eel.is/c++draft/range.factories), which are utilities to create a view.
 *
 * **Sink-only views:** The opposite of a *source-only view*. It can only be placed at the end of a pipe, i.e.
 * it operates on views, but does not actually return a range (has no `rrng_t` column in the second table).
 * The C++20 standard doesn't have this concept.
 *
 * **Deep views:** Some views are declared as "deeps views". This means, that in case they are given a range-of-range
 * as input (as opposed to just a range), they will apply their transformation on the innermost range (instead of
 * the outermost range which would be default). Most alphabet-based transformations are defined as deep, but
 * you can use seqan3::views::deep to make any view (adaptor) deep. See seqan3::views::deep for more details.
 * The C++20 standard doesn't have this concept.
 *
 * **For all views the following are documented:**
 *
 * | Concepts and traits              | `urng_t` (underlying range type) | `rrng_t` (returned range type)                     |
 * |----------------------------------|----------------------------------|----------------------------------------------------|
 * | std::ranges::input_range         | [required] <i>(usually)</i>      | [preserved\|lost\|guaranteed] (usually preserved)  |
 * | std::ranges::forward_range       | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::bidirectional_range | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::random_access_range | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::contiguous_range    | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed] (usually lost)       |
 * |                                  |                                  |                                                    |
 * | std::ranges::viewable_range      | [required] <i>(usually)</i>      | [preserved\|lost\|guaranteed] (usually guaranteed) |
 * | std::ranges::view                | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed] (usually guaranteed) |
 * | std::ranges::sized_range         | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::common_range        | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::output_range        | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | seqan3::const_iterable_range     | [required] <i>(or not)</i>       | [preserved\|lost]                                  |
 * |                                  |                                  |                                                    |
 * | std::ranges::range_reference_t   | optionally a type or concept     | optionally a type or concept                       |
 *
 * **Underlying range requirements:** All view adaptors that are not *source-only* make certain assumptions about their
 * underlying range.
 * The most basic assumption is that the range satisfies `std::ranges::input_range`, but many have stronger
 * requirements, e.g. `std::ranges::random_access_range`. The concepts in the first block all build up on
 * each other, i.e. requiring one implies requiring those above; the other concepts are mostly independent of each
 * other. Most views also require that the underlying range satisfy std::ranges::viewable_range which means they don't
 * accept temporary range objects other than views (because they are cheap to copy). *Note that these being*
 * requirements *means that they are the minimal set of properties assumed. Views may very well make use of stronger
 * properties if available.*
 *
 * **Return range guarantees:** All view adaptors that are not *sink-only* return a range that meets at least
 * `std::ranges::input_range` and also `std::ranges::view` (and conversely also `std::ranges::viewable_range`,
 * because all views are viewable). Most views also preserve stronger
 * properties, e.g. `std::ranges::random_access_range`, but this depends on the view. Some views may add
 * properties not present on the input range.
 *   * *preserved* in this context means that the returned range satisfies this concept if it is also satisfied by the
 * underlying range.
 *   * *lost* means that this concept is never satisfied by the returned range, independent of whether the underlying
 * range supports it.
 *   * *guaranteed* means that this concept is always satisfied by the returned range, independent of whether the underlying
 * range supports it.
 *
 * **Underlying range's reference type:** The reference type is the type the elements of the underlying range are
 * accessed by
 * (since dereferencing an iterator or calling operator[] returns the reference type). The reference type may or may
 * not actually contain a `&` (see below). For many SeqAn specific views additional concept requirements are defined
 * for the input range's reference type, e.g. seqan3::views::complement can only operate on ranges whose elements are
 * nucleotides (meet seqan3::nucleotide_alphabet). In some case the type may even be a specific type or the result
 * of a type trait.
 *
 * **Returned range's reference type:** Conversely certain views make guarantees on the concepts satisfied by the
 * return range's reference type or even always have a fixed type, e.g. seqan3::views::complement operates on
 * nucleotides and of course also returns nucleotides and "std::ranges::range_reference_t<urng_t>" would imply that
 * the reference type is the same. However, and this is important to note, the reference type
 * of seqan3::views::complement has any actual `&` removed from the underlying ranges' reference type (if originally present),
 * this goes hand-in-hand with std::ranges::output_range being lost → original elements cannot be written to through
 * this view.
 * This is because *new elements* are being generated. Other views like `std::views::reverse` also preserve the
 * `&` (if originally present), because the elements in the return view still point to the elements in the original
 * range (just in different order). This has the effect that through some combinations of views you can modify the
 * elements in the original range (if all views in the pipe preserve std::ranges::output_range), but through others
 * you can't.
 *
 * \sa https://en.cppreference.com/w/cpp/ranges
 */

/*!
 * \namespace seqan3::views
 * \brief The SeqAn namespace for views.
 *
 * Since views often have name clashes with regular functions and ranges they are implemented in the sub
 * namespace `view`.
 *
 * See the \link views views submodule \endlink of the range module for more details.
 */

#pragma once

#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/deep.hpp>
#include <seqan3/utility/views/elements.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>
#include <seqan3/utility/views/interleave.hpp>
#include <seqan3/utility/views/join_with.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>
#include <seqan3/utility/views/repeat.hpp>
#include <seqan3/utility/views/repeat_n.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/type_reduce.hpp>
#include <seqan3/utility/views/zip.hpp>
