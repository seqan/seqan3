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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::join_ra.
 */

#pragma once

#include <range/v3/view/join.hpp>

#include <seqan3/core/detail/enum_binary_ops.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/range/concept.hpp>

// --------------------------------------------------------------------------
// ENUM join_flags
// --------------------------------------------------------------------------

namespace seqan3
{

//!\brief Flags to be supplied to seqan3::view::join as template parameters.
//!\ingroup view
enum class view_join_flags : uint8_t
{
    DEFAULT = 0, //!<\brief The default flags.
    SPARSE  = 1, //!<\brief Optimise the data structure for few, long sequences.
//!\cond DEV
    LAZY    = 2  //!<\brief Optimise the data structure for cases that only access the beginning.
//!\endcond
};

/*!\cond DEV
 * \brief Add binary operators to the strictly typed enum seqan3::view_join_flags.
 * \ingroup view
 */
template <>
constexpr bool add_enum_binary_ops<view_join_flags> = true;
//!\endcond

} // namespace seqan3

namespace seqan3::detail
{

// --------------------------------------------------------------------------
// FORWARDS
// --------------------------------------------------------------------------

// unconstrained forward declaration of the view
template <typename irng_t, seqan3::view_join_flags flags>
class view_join_ra;

// unconstrained forward declaration of the iterator
template <typename view_join_ra_t>
class view_join_ra_iterator;

// unconstrained forward declaration of the sentinel
template <typename view_join_ra_t>
class view_join_ra_sentinel;

} // namespace seqan3::detail

#include <seqan3/range/view/join_detail_eager.hpp>
#include <seqan3/range/view/join_detail_lazy.hpp>

namespace ranges::v3
{

template <typename irng_t, seqan3::view_join_flags flags>
struct enable_view<seqan3::detail::view_join_ra<irng_t, flags>> : std::true_type
{};

} // namespace ranges::v3

namespace seqan3::detail
{

// --------------------------------------------------------------------------
// class join_fn (the generator type of view_join_ra)
// --------------------------------------------------------------------------

/*!\brief The type of seqan3::view::join_ra, a generator of seqan3::detail::view_join_ra.
 * \ingroup view
 */
template <seqan3::view_join_flags flags>
struct join_fn
{
    /*!\brief Forwards to ranges::view::join.
     * \tparam irng_t The type of the range being processed. Required to be a seqan3::input_range_concept of
     * seqan3::input_range_concept.
     * \param irange The range being processed.
     */
    template <typename irng_t>
    //!\cond
        requires input_range_concept<irng_t> &&
                 input_range_concept<ranges::range_reference_t<irng_t>>
    //!\endcond
    auto //view_join_ra<ranges::range_reference_t<ranges::range_reference_t<std::decay_t<irng_t>>>>
    operator()(irng_t && irange) const
    {
        return irange | ranges::view::join;
    }

    /*!\brief Special overload for seqan3::concatenated_sequences that just calls
     * seqan3::concatenated_sequences::concat().
     * \tparam value_t Template parameters of seqan3::concatenated_sequences.
     * \tparam delimiter_t Template parameters of seqan3::concatenated_sequences.
     * \param irange The range being processed.
     */
    template <typename value_t, typename delimiter_t>
    auto operator()(concatenated_sequences<value_t, delimiter_t> & irange) const
    {
        return irange.concat();
    }

    /*!\brief Special overload for seqan3::concatenated_sequences that just calls
     * seqan3::concatenated_sequences::concat().
     * \tparam value_t Template parameters of seqan3::concatenated_sequences.
     * \tparam delimiter_t Template parameters of seqan3::concatenated_sequences.
     * \param irange The range being processed.
     */
    template <typename value_t, typename delimiter_t>
    auto operator()(concatenated_sequences<value_t, delimiter_t> const & irange) const
    {
        return irange.concat();
    }

    /*!\brief Overload if input range is not sized. Forwards to the constructor of seqan3::detail::view_join_ra,
     * but add the LAZY flag.
     * \tparam irng_t The type of the range being processed. See seqan3::view::join_ra for requirements.
     * \param irange The range being processed.
     */
    template <typename irng_t>
    //!\cond
        requires random_access_range_concept<irng_t> &&
                 random_access_range_concept<ranges::range_reference_t<irng_t>> &&
                 sized_range_concept<ranges::range_reference_t<irng_t>>
    //!\endcond
    auto //view_join_ra<ranges::range_reference_t<ranges::range_reference_t<std::decay_t<irng_t>>>>
    operator()(irng_t && irange) const
    {
        // if input range is not sized, we must take lazy
        return view_join_ra<std::remove_reference_t<irng_t>,
                            flags | view_join_flags::LAZY>{std::forward<irng_t>(irange)};
    }

    /*!\brief Forwards to the constructor of seqan3::detail::view_join_ra.
     * \tparam irng_t The type of the range being processed. See seqan3::view::join_ra for requirements.
     * \param irange The range being processed.
     */
    template <typename irng_t>
    //!\cond
        requires random_access_range_concept<irng_t> &&
                 sized_range_concept<irng_t> &&
                 random_access_range_concept<ranges::range_reference_t<irng_t>> &&
                 sized_range_concept<ranges::range_reference_t<irng_t>>
    //!\endcond
    auto //view_join_ra<ranges::range_reference_t<ranges::range_reference_t<std::decay_t<irng_t>>>>
    operator()(irng_t && irange) const
    {
        return view_join_ra<std::remove_reference_t<irng_t>, flags>{std::forward<irng_t>(irange)};
    }

    /*!\brief Pipe operator that enables view-typical use of pipe notation.
     * \tparam irng_t The type of the range being processed. Required to be a seqan3::input_range_concept of
     * seqan3::input_range_concept.
     * \param irange The range being processed as left argument of the pipe.
     * \param join_ra_ Instance of seqan3::detail::join_fn that triggers this overload.
     */
    template <typename irng_t>
    //!\cond
        requires input_range_concept<irng_t> &&
                 input_range_concept<ranges::range_reference_t<irng_t>>
    //!\endcond
    friend auto//view_join_ra<ranges::range_reference_t<ranges::range_reference_t<std::decay_t<irng_t>>>>
    operator|(irng_t && irange, join_fn const & join_ra_)
    {
        return join_ra_(std::forward<irng_t>(irange));
    }
};

} // namespace seqan3::detail

namespace seqan3::view
{

/*!\brief For a range of ranges, return the flattened range.
 * \tparam irng_t The type of the range being processed. See below for requirements.
 * [template parameter is omitted in pipe notation]
 * \tparam flags Flags to fine-tune the view's behaviour, see seqan3::view_join_flags.
 * \param irange The range being processed. [parameter is omitted in pipe notation]
 * \returns A flattened range. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * Similar to ranges::view::join, but with additional features.
 *
 * \cond !DEV
 * |  case |  Description                                         | Effect                                           |
 * |-------|------------------------------------------------------|--------------------------------------------------|
 * |  1.   | *if none of the below*                               | same as ranges::view::join                       |
 * |  2.   | input is seqan3::concatenated_sequences              | return seqan3::concatenated_sequences::concat()  |
 * |  3.   | both dimensions of input have random access and size | preserve random access and size in return type   |
 * |  3.s  | 3. and seqan3::view_join_flag::SPARSE was set        | faster access if sub-ranges are long             |
 * \endcond
 *
 * \cond DEV
 * |  case |  Description                                         | Effect                                           |
 * |-------|------------------------------------------------------|--------------------------------------------------|
 * |  1.   | *if none of the below*                               | same as ranges::view::join                       |
 * |  2.   | input is seqan3::concatenated_sequences              | return seqan3::concatenated_sequences::concat()  |
 * |  3.   | both dimensions of input have random access and size | preserve random access and size in return type   |
 * |  3.s  | 3. and seqan3::view_join_flag::SPARSE was set        | lower space consumption, faster access if seqs   |
 * |  3.l  | 3. and seqan3::view_join_flag::LAZY was set          | faster access if sub-ranges are long             |
 * \endcond
 *
 * Does not provide an interface for inserting delimiters like ranges::view::join.
 *
 * \attention
 *
 * * If you need the best-possible performance, use seqan3::concatenated_sequences for your input data.
 * * If the average size of your sub-ranges is >= 1000, use the sparse specialisation.
 *
 * \par Example
 *
 * Case 1:
 * ```cpp
 * // input is only a forward_range
 * std::forward_list<dna5_vector>      vec{"AAAAA"_dna5, "CCCC"_dna5, "GGG"_dna5, "TT"_dna5};
 *
 * auto v = vec | view::join<>;   // this is the same as "auto v = vec | ranges::view::join;"
 *                                // any template parameters are irrelevant
 *
 * // assert(v[6] == dna5::C);       // not supported
 * // assert(ranges::size(v) == 14); // not supported
 * ```
 *
 * Case 2:
 * ```cpp
 * // input is seqan3::concatenated_sequences
 * concatenated_sequences<dna5_vector> vec{"AAAAA"_dna5, "CCCC"_dna5, "GGG"_dna5, "TT"_dna5};
 *
 * auto v = vec | view::join<>;   // this is the same as "auto v = vec.concat();"
 *                                // any template parameters are irrelevant
 *
 * assert(v[6] == dna5::C);       // supported and is the fastest possible
 * assert(ranges::size(v) == 14); // supported
 * ```
 *
 * Case 3:
 * ```cpp
 * // input is a sized random access range, e.g. a vector
 * std::vector<dna5_vector>            vec{"AAAAA"_dna5, "CCCC"_dna5, "GGG"_dna5, "TT"_dna5};
 *
 * auto v = vec | view::join<>;   // this is the same as "auto v = vec | view::join<view_join_flags::DEFAULT>;"
 *
 * // use this if average sub-range length >= 1000:
 * // auto v = vec | view::join<view_join_flags::SPARSE>;
 *
 * assert(v[6] == dna5::C);       // supported, but not as fast as on seqan3::concatenated_sequences
 * assert(ranges::size(v) == 14); // supported
 * ```
 *
 * \par View properties
 *
 * * The input properties are **requirements** on the range input type.
 * * The return properties are **guarantees** given on the range return type.
 * * For descriptions of the view properties see \ref view.
 *
 * \cond !DEV
 * | case |                     | `irng_t` (range input type)    | `rrng_t` (range return type)                              |
 * |------|---------------------|--------------------------------|-----------------------------------------------------------|
 * | 1.   | range               | seqan3::input_range_concept    | seqan3::view_concept + seqan3::input_range_concept (all other are lost)       |
 * |      | `range_reference_t` | seqan3::input_range_concept    | `range_reference_t<range_reference_t<irng_t>>`            |
 * | 2.   | range               | seqan3::concatenated_sequences | seqan3::view_concept + seqan3::random_access_range_concept + seqan3::sized_range_concept |
 * |      | `range_reference_t` | seqan3::random_access_range_concept + seqan3::sized_range_concept   | `range_reference_t<range_reference_t<irng_t>>`      |
 * | 3.*  | range               | seqan3::random_access_range_concept + seqan3::sized_range_concept  | seqan3::view_concept + seqan3::random_access_range_concept + seqan3::sized_range_concept |
 * |      | `range_reference_t` | seqan3::random_access_range_concept + seqan3::sized_range_concept  | `range_reference_t<range_reference_t<irng_t>>`      |
 *
 * | case  | Allocating view   |
 * |-------|:-----------------:|
 * | 1.    | ☐                 |
 * | 2.    | ☐                 |
 * | 3.*   | ☑                 |
 * \endcond
 *
 * \cond DEV
 * | case |                     | `irng_t` (range input type)    | `rrng_t` (range return type)                              |
 * |------|---------------------|--------------------------------|-----------------------------------------------------------|
 * | 1.   | range               | seqan3::input_range_concept    | seqan3::view_concept + seqan3::input_range_concept (all other are lost)       |
 * |      | `range_reference_t` | seqan3::input_range_concept    | `range_reference_t<range_reference_t<irng_t>>`            |
 * | 2.   | range               | seqan3::concatenated_sequences | seqan3::view_concept + seqan3::random_access_range_concept + seqan3::sized_range_concept |
 * |      | `range_reference_t` | seqan3::random_access_range_concept + seqan3::sized_range_concept   | `range_reference_t<range_reference_t<irng_t>>`      |
 * | 3. not l | range               | seqan3::random_access_range_concept + seqan3::sized_range_concept  | seqan3::view_concept + seqan3::random_access_range_concept + seqan3::sized_range_concept |
 * |      | `range_reference_t` | seqan3::random_access_range_concept + seqan3::sized_range_concept  | `range_reference_t<range_reference_t<irng_t>>`      |
 * | 3.l  | range               | seqan3::random_access_range_concept | seqan3::view_concept + seqan3::random_access_range_concept |
 * |      | `range_reference_t` | seqan3::random_access_range_concept + seqan3::sized_range_concept  | `range_reference_t<range_reference_t<irng_t>>`      |
 *
 * | case  | Allocating view   | Const-Iterable  |
 * |-------|:-----------------:|:---------------:|
 * | 1.    | ☐                 |      ☑          |
 * | 2.    | ☐                 |      ☑          |
 * | 3.    | ☑                 |      ☑          |
 * | 3.s   | ☑                 |      ☑          |
 * | 3.l   | ☑                 |      ☐          |
 * | 3.s+l | ☑                 |      ☐          |
 *\endcond
 *
 * \par Complexity
 *
 * * let \f$m\f$ be the number of elements (sub-ranges) in `irange`
 * * let \f$n\f$ be the total joined length
 *
 * \cond !DEV
 * |  case  | creation   | `[]` (anywhere)                   | `[]` (same sub-range)   | space complexity[bits] |
 * |--------|------------|-----------------------------------|-------------------------|------------------------|
 * |  1.    | \f$O(1)\f$ | --                                | --                      | \f$O(1)\f$             |
 * |  2.    | \f$O(1)\f$ | \f$O(1)\f$                        | \f$O(1)\f$              | \f$O(1)\f$             |
 * |  3.    | \f$O(n)\f$ | \f$O(\log(^n/_m))\f$              | \f$O(1)\f$              | \f$\approx m * (2 + \log(^n/_m))\f$ |
 * |  3.s   | \f$O(m)\f$ | \f$O(\log(m))\f$                  | \f$O(1)\f$              | \f$\approx 64*m\f$     |
 * \endcond
 *
 * \cond DEV
 * |  case  | creation   | `[]` (anywhere)                   | `[]` (same sub-range)   | space complexity[bits] |
 * |--------|------------|-----------------------------------|-------------------------|------------------------|
 * |  1.    | \f$O(1)\f$ | --                                | --                      | \f$O(1)\f$             |
 * |  2.    | \f$O(1)\f$ | \f$O(1)\f$                        | \f$O(1)\f$              | \f$O(1)\f$             |
 * |  3.    | \f$O(n)\f$ | \f$O(\log(^n/_m))\f$              | \f$O(1)\f$              | \f$\approx m * (2 + \log(^n/_m))\f$ |
 * |  3.s   | \f$O(m)\f$ | \f$O(\log(m))\f$                  | \f$O(1)\f$              | \f$\approx 64*m\f$     |
 * |  3.l   | \f$O(1)\f$ | *amortised* \f$O(\log(^n/_m)))\f$ | \f$O(1)\f$              | \f$\approx m * (2 + \log(^n/_m))\f$ |
 * |  3.l+s | \f$O(1)\f$ | *amortised* \f$O(\log(m))\f$      | \f$O(1)\f$              | \f$\approx 64*m\f$     |
 * \endcond
 *
 * As mentioned above, it is recommended to use the sparse version once the average sub-range length (\f$^n/_m\f$)
 * is over 1000.
 */

template <view_join_flags flags = view_join_flags::DEFAULT>
seqan3::detail::join_fn<flags> const join;

} // namespace seqan3::view

#ifndef NDEBUG
#include <seqan3/range/view/concept.hpp>
//!\cond
#define PARAM_TYPES std::vector<std::vector<char>>, seqan3::view_join_flags::DEFAULT
//!\endcond
static_assert(seqan3::random_access_iterator_concept<
                seqan3::detail::view_join_ra_iterator<seqan3::detail::view_join_ra<PARAM_TYPES> const>>);

static_assert(seqan3::input_range_concept<seqan3::detail::view_join_ra<PARAM_TYPES>>);
static_assert(seqan3::random_access_range_concept<seqan3::detail::view_join_ra<PARAM_TYPES>>);
static_assert(seqan3::view_concept<seqan3::detail::view_join_ra<PARAM_TYPES>>);
#undef PARAM_TYPES
#endif
