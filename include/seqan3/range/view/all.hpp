// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link view view submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/interleave.hpp>
#include <seqan3/range/view/istreambuf.hpp>
#include <seqan3/range/view/pairwise_combine.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/rank_to.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/range/view/translate.hpp>
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
 * See the \link range range module \endlink for how views relate to containers and decorators.
 *
 * Most views provided by SeqAn3 are specific to biological operations, like seqan3::view::trim which trims
 * sequences based on the quality or seqan3::view::complement which generates the complement of a nucleotide sequence.
 * But SeqAn3 also provides some general purpose views.
 *
 * ### Namespaces
 *
 *   * [All views from the range-v3 libary](https://ericniebler.github.io/range-v3/index.html#range-views) are available
 * in the namespace `ranges::view`.
 *
 *   * All SeqAn views are available in the namespace `seqan3::view`.
 *
 * ### Example
 *
 * Functional and pipe notations:
 * \snippet test/snippet/range/view/range_view_all.cpp notation
 *
 * Re-transform into a distinct container:
 * \snippet test/snippet/range/view/range_view_all.cpp retransform
 *
 * Composability:
 * \snippet test/snippet/range/view/range_view_all.cpp composability
 *
 * ### Views vs view adaptors
 *
 * When talking about views, two different entities are often conflated:
 *
 *   1. the view (this is the type that is a range and meets std::ranges::View; it is what we refer to with
 * `auto vec_view` above)
 *   2. the view adaptor (this is the functor that returns the actual view based on it's parameters, including the
 * underlying range; in the above example `view::reverse` and `view::complement` are view adaptors)
 *
 * The view adaptor also facilitates the piping behaviour. It is the only entity that is publicly documented and
 * the actual view type (the range type returned by the adaptor) is considered implementation defined.
 * The *properties* of the returned range are however specified and documented as part of the adaptor, see below.
 *
 * <sub>
 * An exception to this rule are views that don't work on an underlying range and can only be
 * placed at the beginning of a pipe of operations; they do not need an adaptor, because their constructor is
 * sufficient. This is not relevant for the documentation, though, we always document `view::foo`, independent of
 * whether `view::foo` is the adaptor that returns the "foo type" or whether `view::foo` is the "foo type".
 * </sub>
 *
 * ### View properties
 *
 * There are three view properties that are documented for a view, **only if that view fulfills them:**
 *
 * **Source-only views:** Most views operate on an underlying range and return a (modified) range, i.e. they can be placed
 * at the beginning, middle or end of a "pipe" of view operations. However, some views are limited to being at
 * the front ("source"), e.g. `std::view::single`, `ranges::view::concat` and `ranges::view::ints`. These views
 * are marked as "source-only" and have no `urng_t` column in the second table.
 *
 * **Sink-only views:** The opposite of a *source-only view*. It can only be placed at the end of a pipe, i.e.
 * it operates on views, but does not actually return a range (has no `rrng_t` column in the second table).
 *
 * **Deep views:** Some views are declared as "deeps views". This means, that in case they are given a range-of-range
 * as input (as opposed to just a range), they will apply their transformation on the innermost range (instead of
 * the outermost range which would be default). Most alphabet-based transformations are defined as deep, but
 * you can use seqan3::view::deep to make any view (adaptor) deep. See seqan3::view::deep for more details.
 *
 * **For all views the following are documented:**
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type) | `rrng_t` (returned range type)                     |
 * |---------------------------------|----------------------------------|----------------------------------------------------|
 * | std::ranges::InputRange         | [required] <i>(usually)</i>      | [preserved\|lost\|guaranteed] (usually preserved)  |
 * | std::ranges::ForwardRange       | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::BidirectionalRange | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::RandomAccessRange  | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::ContiguousRange    | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed] (usually lost)       |
 * |                                 |                                  |                                                    |
 * | std::ranges::ViewableRange      | [required] <i>(usually)</i>      | [preserved\|lost\|guaranteed] (usually guaranteed) |
 * | std::ranges::View               | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed] (usually guaranteed) |
 * | std::ranges::SizedRange         | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::CommonRange        | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | std::ranges::OutputRange        | [required] <i>(or not)</i>       | [preserved\|lost\|guaranteed]                      |
 * | seqan3::ConstIterableRange      | [required] <i>(or not)</i>       | [preserved\|lost]                                  |
 * |                                 |                                  |                                                    |
 * | seqan3::reference_t             | optionally a type or concept     | optionally a type or concept                       |
 *
 * **Underlying range requirements:** All view adaptors that are not *source-only* make certain assumptions about their
 * underlying range.
 * The most basic assumption is that the range satisfies `std::ranges::InputRange`, but many have stronger
 * requirements, e.g. `std::ranges::RandomAccessRange`. The concepts in the first block all build up on
 * each other, i.e. requiring one implies requiring those above; the other concepts are mostly independent of each
 * other. Most views also require that the underlying range satisfy std::ranges::ViewableRange which means they don't
 * accept temporary range objects other than views (because they are cheap to copy). A prominent exception
 * to the latter is view::persist that exists exactly for this purpose. *Note that these being* requirements *means that
 * they are the minimal set of properties assumed. Views may very well make use of stronger properties if available.*
 *
 * **Return range guarantees:** All view adaptors that are not *sink-only* return a range that meets at least
 * `std::ranges::InputRange` and also `std::ranges::View` (and conversely also `std::ranges::ViewableRange`,
 * because all views are viewable). Most views also preserve stronger
 * properties, e.g. `std::ranges::RandomAccessRange`, but this depends on the view. Some views also add
 * properties not present on the input range, e.g. the range returned by `ranges::view::take_exactly` meets
 * `std::ranges::SizedRange`, independent of whether this was met by the input range.
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
 * for the input range's reference type, e.g. seqan3::view::complement can only operate on ranges whose elements are
 * nucleotides (meet seqan3::NucleotideAlphabet_check). In some case the type may even be a specific type or the result
 * of a type trait.
 *
 * **Returned range's reference type:** Conversely certain views make guarantees on the concepts satisfied by the
 * return range's reference type or even always have a fixed type, e.g. seqan3::view::complement operates on
 * nucleotides and of course also returns nucleotides and "seqan3::reference_t<urng_t>" would imply that
 * the reference type is the same. However, and this is important to note, the reference type
 * of seqan3::view::complement has any actual `&` removed from the underlying ranges' reference type (if originally present),
 * this goes hand-in-hand with std::ranges::OutputRange being lost → original elements cannot be written to through
 * this view.
 * This is because *new elements* are being generated. Other views like `view::reverse` also preserve the
 * `&` (if originally present), because the elements in the return view still point to the elements in the original
 * range (just in different order). This has the effect that through some combinations of views you can modify the
 * elements in the original range (if all views in the pipe preserve std::ranges::OutputRange), but through others
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
