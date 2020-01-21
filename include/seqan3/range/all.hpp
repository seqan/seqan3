// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link range range module \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/container/all.hpp>
#include <seqan3/range/decorator/all.hpp>
#include <seqan3/range/views/all.hpp>

/*!\defgroup range Range
 * \brief The range module provides general purpose containers, decorators and views.
 *
 * ### Introduction
 *
 * *Ranges* are an abstraction of "a collection of items", or "something iterable". The most basic definition
 * requires only the existence of `begin()` and `end()` on the range.
 *
 * There are different ways to classify ranges, one way is through the capabilities of its default iterators.
 * This is resembled by the range concepts defined in this module. Another way to classify ranges is by their storage
 * behaviour, i.e. whether they own the data that is accessible through them. See below for more details.
 *
 * Ranges are found throughout the SeqAn library, this module provides general-purpose ranges that are not specific
 * to another module or biological function.
 *
 * ### Iterator capabilities
 *
 * All ranges in SeqAn are either \link std::ranges::input_range input ranges \endlink (they can be read from) or
 * \link std::ranges::output_range output ranges \endlink (they can be written to) or both. E.g. an
 * `std::vector<int>` is both, but a `std::vector<int> const` would only be an input range.
 *
 * \link std::ranges::input_range Input ranges \endlink have different *strengths* that are realised through more
 * refined concepts:
 *
 * std::ranges::input_range < std::ranges::forward_range < std::ranges::bidirectional_range <
 * std::ranges::random_access_range < std::ranges::contiguous_range
 *
 * (Click on the respective concepts to learn the exact definitions)
 *
 * Independent of input or output, a range can also be \link std::ranges::sized_range sized \endlink and/or
 * \link std::ranges::common_range  common \endlink.
 *
 * ### Storage behaviour
 *
 * **Containers** are the ranges most well known, they own their elements. SeqAn makes use of standard STL containers
 * like `std::vector`, but also implements some custom containers. See the \link container container submodule \endlink
 * for more details.
 *
 * **Decorators** are ranges that are always defined on another range and decorate/annotate the underlying range
 * with additional information. They do not own the underlying range, but can contain member data of their own.
 * See the \link decorator decorator submodule \endlink for more details.
 *
 * **Views** are ranges that are usually defined on another range and transform the underlying range
 * via some algorithm or operation, however, some views are stand-alone, i.e. they are just an algorithm that
 * produces elements. Views do not own any data beyond their algorithm and possible parameters to it and they
 * can always be copied in constant time. The algorithm is required to be lazy-evaluated so it is feasible to
 * combine multiple views. See the \link views views submodule \endlink for more details.
 *
 * If you are confused about *decorators* vs *views*, think of decorators as "underlying range + data" and
 * views as "underlying range + algorithm".
 *
 * The storage behaviour is orthogonal to the range concepts defined by the iterators mentioned above, i.e. you
 * can have a container that satisfies std::ranges::random_access_range (e.g. `std::vector` does, but `std::list`
 * does not) and you can have views or decorators that do so or don't. For some combinations of iterator capabilities
 * and storage behaviour there are extra concept definitions, e.g. seqan3::random_access_container.
 *
 * \attention
 *
 * There are ranges in SeqAn that fit neither of these storage categories, e.g. all the files are
 * \link std::ranges::input_range input ranges \endlink (if they are input files) and
 * \link std::ranges::output_range output ranges \endlink (if they are output files), but they are neither
 * containers, decorators nor views.
 *
 * \sa range.hpp
 * \sa https://ericniebler.github.io/range-v3/index.html
 */
