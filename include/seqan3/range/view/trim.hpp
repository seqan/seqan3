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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::trim.
 */

#pragma once

#include <range/v3/view/take_while.hpp>

#include <seqan3/alphabet/quality/quality_composition.hpp>
#include <seqan3/std/concept/range.hpp>
#include <seqan3/range/view/deep.hpp>

namespace seqan3::detail
{

/*!\brief The underlying type of seqan3::view::trim.
 * \ingroup view
 *
 * Under the hood this delegates to ranges::view::take_while.
 */
struct trim_fn
{
    /*!\brief Trim based on minimum phred score.
     * \tparam irng_t The type of the range being processed. See seqan3::view::trim for requirements.
     * \param irange The range being processed.
     * \param threshold The minimum quality as a phred score [integral type].
     */
    template <typename irng_t>
    auto operator()(irng_t && irange,
                    underlying_phred_t<ranges::value_type_t<std::decay_t<irng_t>>> const threshold) const
    //!\cond
        requires input_range_concept<irng_t> && quality_concept<ranges::value_type_t<std::decay_t<irng_t>>>
    //!\endcond
    {
        return ranges::view::take_while(std::forward<irng_t>(irange), [threshold] (auto && value)
        {
            return to_phred(std::forward<decltype(value)>(value)) >= threshold;
        });
    }

    /*!\brief Trim based on value_type.
     * \tparam irng_t The type of the range being processed. See seqan3::view::trim for requirements.
     * \param irange The range being processed.
     * \param threshold The minimum quality given by a value of the ranges type.
     */
    template <typename irng_t>
    auto operator()(irng_t && irange,
                    std::decay_t<ranges::value_type_t<std::decay_t<irng_t>>> const threshold) const
    //!\cond
        requires input_range_concept<irng_t> && quality_concept<ranges::value_type_t<std::decay_t<irng_t>>>
    //!\endcond
    {
        return (*this)(std::forward<irng_t>(irange), to_phred(threshold));
    }

    /*!\brief A functor that behaves like a named version of std::bind around seqan3::detail::trim_fn::operator().
     * \tparam threshold_t Must be an integral type or satisfy the seqan3::quality_concept.
     *
     * \details
     *
     * * The single-parameter-operator() of trim_fn would normally return a bound two-parameter-operator() via
     * std::bind.
     * * trim_fn's friend operator| would take this as second argument via decltype(std::bind(...)) so that it
     * is specific to this overload.
     * * However, since the type of the threshold parameter is a template argument, we need function template
     * overloading. This doesn't work with decltype() and the actual return type of std::bind is implementation
     * defined.
     * * So we have this helper class template which behaves exactly the same, but has a distinct (named) type.
     *
     * \attention You should never instantiate this manually.
     */
    template <typename threshold_t>
    struct delegate
    {
        //!\brief The intermediately stored threshold.
        threshold_t const threshold;
        //!\brief Reference to the parent
        trim_fn const & parent;

        //!\brief The operator() that only takes the range as argument and forwards to the two-parameter operator().
        template <typename irng_t>
        auto operator()(irng_t && irange) const
        {
            return parent(std::forward<irng_t>(irange), threshold);
        }
    };

    /*!\brief Range-less interface for use with the pipe notation.
     * \tparam threshold_t Must be an integral type or satisfy the seqan3::quality_concept.
     * \param threshold The minimum quality given by a value of the range's value type.
     *
     * \details
     *
     * Binds to one of the other interfaces and forwards the threshold.
     */
    template <typename threshold_t>
    delegate<threshold_t> operator()(threshold_t const threshold) const
    //!\cond
        requires std::is_integral_v<std::decay_t<threshold_t>> || quality_concept<std::decay_t<threshold_t>>
    //!\endcond
    {
        return delegate<threshold_t>{threshold, *this};
        // this doesn't work here, see seqan3::detail::trim_fn::delegate.
        //return std::bind(trim_fn(), std::placeholders::_1, threshold);
    }

    /*!\brief Pipe operator that enables view-typical use of pipe notation.
     * \tparam irng_t The type of the range being processed. See seqan3::view::trim for requirements.
     * \tparam threshold_t Must be of `irng_t`'s `value_type` or that `value_type`'s seqan3::underlying_phred_t.
     * \param irange The range being processed as left argument of the pipe.
     * \param bound_view The result of the single-argument operator() (interface with bound threshold parameter).
     */
    template <typename irng_t,
              typename threshold_t>
    //!\cond
        requires input_range_concept<irng_t> && quality_concept<ranges::value_type_t<std::decay_t<irng_t>>> &&
                 (std::is_same_v<std::decay_t<threshold_t>,
                                 std::decay_t<ranges::value_type_t<std::decay_t<irng_t>>>> ||
                  std::is_convertible_v<std::decay_t<threshold_t>,
                                        underlying_phred_t<std::decay_t<ranges::value_type_t<std::decay_t<irng_t>>>>>)
    //!\endcond
    friend auto operator|(irng_t && irange,
                          seqan3::detail::trim_fn::delegate<threshold_t> const & bound_view)
    {
        return bound_view(std::forward<irng_t>(irange));
    }
};

} // namespace seqan3::detail

namespace seqan3::view
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               A view that does quality-threshold trimming on a range of seqan3::quality_concept.
 * \tparam urng_t       The type of the range being processed. See below for requirements.
 * \tparam threshold_t  Either seqan3::value_type_t<urng_t> or
 *                      seqan3::underlying_phred_t<seqan3::value_type_t<urng_t>>.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] threshold The minimum quality.
 * \returns             A trimmed range. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view can be used to do easy quality based trimming of sequences.
 *
 * ### View properties
 *
 * This view is a **deep view:** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | range concepts and reference_t      | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |-------------------------------------|:-------------------------------------:|:-------------------------------:|
 * | seqan3::input_range_concept         | *required*                            | *preserved*                     |
 * | seqan3::forward_range_concept       |                                       | *preserved*                     |
 * | seqan3::bidirectional_range_concept |                                       | *preserved*                     |
 * | seqan3::random_access_range_concept |                                       | *preserved*                     |
 * |                                     |                                       |                                 |
 * | seqan3::view_concept                |                                       | *guaranteed*                    |
 * | seqan3::sized_range_concept         |                                       | *lost*                          |
 * | seqan3::bounded_range_concept       |                                       | *lost*                          |
 * | seqan3::output_range_concept        |                                       | *preserved*                     |
 * | seqan3::const_iterable_concept      |                                       | *preserved*                     |
 * |                                     |                                       |                                 |
 * | seqan3::reference_t                 | seqan3::quality_concept               | seqan3::reference_t<urng_t>     |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * \par Example
 *
 * Operating on a range of seqan3::illumina18:
 * ```cpp
 * std::vector<illumina18> vec{illumina18{40}, illumina18{40}, illumina18{30}, illumina18{20}, illumina18{10}};
 *
 * // trim by phred_value
 * auto v1 = vec | view::trim(20u);                        // == ['I','I','?','5']
 *
 * // trim by quality character
 * auto v2 = vec | view::trim(illumina18{40});             // == ['I','I']
 *
 * // function syntax
 * auto v3 = view::trim(vec, 20u);                         // == ['I','I','?','5']
 *
 * // combinability
 * std::string v4 = view::trim(vec, 20u) | view::to_char;  // == "II?5"
 * ```
 *
 * Or operating on a range of seqan3::dna5q:
 * ```cpp
 * std::vector<dna5q> vec{{dna5::A, illumina18{40}}, {dna5::G, illumina18{40}}, {dna5::G, illumina18{30}},
 *                        {dna5::A, illumina18{20}}, {dna5::T, illumina18{10}}};
 * std::vector<dna5q> cmp{{dna5::A, illumina18{40}}, {dna5::G, illumina18{40}}, {dna5::G, illumina18{30}},
 *                        {dna5::A, illumina18{20}}};
 *
 * // trim by phred_value
 * auto v1 = vec | view::trim(20u);
 * assert(std::vector<dna5q>(v1) == cmp);
 *
 * // trim by quality character; in this case the nucleotide part of the character is irrelevant
 * auto v2 = vec | view::trim(dna5q{dna5::C, illumina18{20}});
 * assert(std::vector<dna5q>(v2) == cmp);
 *
 * // combinability
 * std::string v4 = view::trim(vec, 20u) | view::to_char;
 * EXPECT_EQ("AGGA", v4);
 * ```
 * \hideinitializer
 */

inline constexpr auto trim = deep{seqan3::detail::trim_fn{}};

//!\}

} // namespace seqan3::view
