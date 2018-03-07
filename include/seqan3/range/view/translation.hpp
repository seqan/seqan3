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
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::view::translate, seqan3::view::translate_single and seqan3::view::translate_frames.
 */

#pragma once

#include <stdexcept>
#include <vector>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation.hpp>
#include <seqan3/alphabet/aminoacid/translation_details.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/complement.hpp>

#include <range/v3/view/chunk.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/take_while.hpp>
#include <range/v3/view/transform.hpp>

namespace seqan3
{
/*!\brief Specialization values for translation frames.
 * \details The numeric values of the enums correspond to the number of frames that are represented.
 */
enum class translation_frames
{
    SINGLE_FRAME = 1,
    WITH_REVERSE_COMPLEMENT = 2,
    WITH_FRAME_SHIFTS = 3,
    SIX_FRAME = 6
};
}

namespace seqan3::detail
{
/*!\brief Lambda for translation of one nucleotide triplet into a single aminoacid.
 * \details Used by view::translate_single.
 */
auto translate = [](auto range) { return translate_triplet(range); };

/*!\brief The underlying type of seqan3::view::translate_single.
 * \ingroup view
 */
struct translate_single_fn
{
    /*!\brief Single frame translation from nucleotide to aminoacid.
     * \tparam irng_t The type of the range being processed.
     * \param irange The range being processed.
     */
    template <typename irng_t>
    //!\cond
    requires input_range_concept<irng_t> &&
      nucleotide_concept<std::decay_t<ranges::range_reference_t<std::decay_t<irng_t>>>>
      //!\endcond
      auto operator()(irng_t && irange) const
    {
        if constexpr (sized_range_concept<irng_t>) {
            return irange | ranges::view::take_exactly(ranges::size(irange) / 3 * 3) | ranges::view::chunk(3) |
                   ranges::view::transform(translate);
        }
        else
        {
            return irange | ranges::view::chunk(3) |
                   ranges::view::take_while([](auto chnk) { return ranges::size(chnk) == 3; }) |
                   ranges::view::transform(translate);
        }
    }

    /*!\brief Pipe operator that enables view-typical use of pipe notation.
    * \tparam irng_t The type of the range being processed.
    * \param irange The range being processed as left argument of the pipe.
    * \param translate_single_fn_v The result of the single-argument operator().
    */
    template <typename irng_t>
    friend auto operator|(irng_t && irange, seqan3::detail::translate_single_fn const & translate_single_fn_v)
    {
        return translate_single_fn_v(std::forward<irng_t>(irange));
    }
};

} // namespace seqan3::detail

namespace seqan3::view
{
/*!\brief A view performs single frame translation of nucleotide into aminoacid alphabet.
 * \tparam irng_t The type of the range being processed.
 * \param irange The range being processed.
 * \returns A range containing aminoacid sequence. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view can be used to translate nucleotide sequences into single forward aminoacid sequences.

 * \par View properties
 *
 * |                     | `irng_t` (range input type)           | `rrng_t` (range return type)                              |
 * |---------------------|---------------------------------------|-----------------------------------------------------------|
 * | range               | seqan3::input_range_concept           | seqan3::view_concept + all range concepts met by `irng_t` |
 * | `range_reference_t` | seqan3::nucleotide_concept            |                                                           |
 *
 * * `irng_t` is the type of the range modified by this view (input).
 * * `rrng_type` is the type of the range returned by this view.
 * * for more details, see \ref view.
 *
 * \par Example
 *
 * Operating on a range of seqan3::dna5:
 * ```cpp
 * dna5_vector vec{"ACGTACGTACGTA"_dna5};
 *
 * // single frame translation
 * auto v1 = vec | view::translate_single;                     // == [T,Y,V,R]
 *
 * // function syntax
 * auto v2 = view::translate_single(vec);                      // == [T,Y,V,R]
 *
 * // combinability
 * auto v3 = vec | view::complement | view::translate_single;  // == [C,M,H,A]
 * ```
 */

seqan3::detail::translate_single_fn const translate_single;

} // namespace seqan3::view

namespace seqan3::detail
{
/*!\brief The underlying type of seqan3::view::translate_frames.
 * \ingroup view
 */
struct translate_frames_fn
{
    /*!\brief Single, reverse, forward-frame or six-frame translation from nucleotide to aminoacid.
     * \tparam irng_t The type of the range being processed.
     * \param irange The range being processed.
     * \param tf Translation frames (translation_frames::SINGLE_FRAME, translation_frames::WITH_REVERSE_COMPLEMENT,
     * translation_frames::WITH_FRAME_SHIFTS, translation_frames::SIX_FRAME).
     */
    template <typename irng_t>
    //!\cond
    requires forward_range_concept<irng_t> &&
      nucleotide_concept<std::decay_t<ranges::range_reference_t<std::decay_t<irng_t>>>>
      //!\endcond
      auto operator()(irng_t && irange, translation_frames const & tf) const
    {
        std::vector<ranges::any_view<aa27, ranges::category::random_access>> frames;
        frames.resize(static_cast<uint8_t>(tf));

        switch (tf)
        {
            case translation_frames::WITH_REVERSE_COMPLEMENT:
                frames[1] = irange | ranges::view::reverse | view::complement | view::translate_single;
                [[fallthrough]];
            case translation_frames::SINGLE_FRAME: frames[0] = irange | view::translate_single; break;
            case translation_frames::SIX_FRAME:
                for (unsigned i : { 0, 1, 2 })
                    frames[i + 3] = irange | ranges::view::reverse | view::complement | ranges::view::drop(i) |
                                    view::translate_single;
                [[fallthrough]];
            case translation_frames::WITH_FRAME_SHIFTS:
                for (unsigned i : { 0, 1, 2 }) frames[i] = irange | ranges::view::drop(i) | view::translate_single;
                break;
            default:
                throw std::invalid_argument("Invalid number of frames. Choose from SINGLE_FRAME, "
                                            "WITH_REVERSE_COMPLEMENT, WITH_FRAME_SHIFTS, SIX_FRAME.");
                break;
        }

        ranges::any_view<ranges::any_view<aa27, ranges::category::random_access>, ranges::category::random_access>
          combined;
        switch (tf)
        {
            case translation_frames::SINGLE_FRAME: combined = ranges::view::single(frames[0]); break;
            case translation_frames::WITH_REVERSE_COMPLEMENT:
                combined = ranges::view::concat(ranges::view::single(frames[0]), ranges::view::single(frames[1]));
                break;
            case translation_frames::WITH_FRAME_SHIFTS:
                combined = ranges::view::concat(
                  ranges::view::single(frames[0]), ranges::view::single(frames[1]), ranges::view::single(frames[2]));
                break;
            case translation_frames::SIX_FRAME:
                combined = ranges::view::concat(ranges::view::single(frames[0]),
                                                ranges::view::single(frames[1]),
                                                ranges::view::single(frames[2]),
                                                ranges::view::single(frames[3]),
                                                ranges::view::single(frames[4]),
                                                ranges::view::single(frames[5]));
                break;
            default: break;
        }
        return combined;
    }

    /*!\brief Range-less interface for use with the pipe notation.
     * \tparam size_type The type of the input number of frames. Must be an integral type.
     * \param tf Translation frames (translation_frames::SINGLE_FRAME, translation_frames::WITH_REVERSE_COMPLEMENT,
     * translation_frames::WITH_FRAME_SHIFTS, translation_frames::SIX_FRAME).
     *
     * \details
     * Binds to translate_frames_fn and forwards the number of frames.
     */
    auto operator()(translation_frames const tf) const
    {
        return std::bind(translate_frames_fn(), std::placeholders::_1, tf);
    }

    /*!\brief Pipe operator that enables view-typical use of pipe notation.
     * \tparam irng_t The type of the range being processed.
     * \param irange The range being processed as left argument of the pipe.
     * \param bound_view The result of the single-argument operator() (interface with bound tf parameter).
     */
    template <random_access_range_concept irng_t>
    friend auto operator|(irng_t && irange,
                          decltype(std::bind(translate_frames_fn(),
                                             std::placeholders::_1,
                                             translation_frames::SIX_FRAME)) const & bound_view)
    {
        return bound_view(std::forward<irng_t>(irange));
    }
};

} // namespace seqan3::detail

namespace seqan3::view
{
/*!\brief A view that translates nucleotide into aminoacid alphabet with 1, 2, 3 or 6 frames.
 * \tparam irng_t The type of the range being processed.
 * \param irange The range being processed.
 * \param tf Translation frames (translation_frames::SINGLE_FRAME, translation_frames::WITH_REVERSE_COMPLEMENT,
 * translation_frames::WITH_FRAME_SHIFTS, translation_frames::SIX_FRAME).
 * \returns A range of ranges containing frames with aminoacid sequence. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view can be used to translate nucleotide sequences into aminoacid sequences:
 * translation_frames::SINGLE_FRAME             - forward frame (Note: The output will be a range of ranges. If a single range
 *                                                is wanted or needed, please use view::translate_single)
 * translation_frames::WITH_REVERSE_COMPLEMENT  - forward/ reverse frames
 * translation_frames::WITH_FRAME_SHIFTS        - all forward frames
 * translation_frames::SIX_FRAME                - all forward and reverse frames
 *
 * \par View properties
 *
 * |                     | `irng_t` (range input type)           | `rrng_t` (range return type)                                |
 * |---------------------|---------------------------------------|-------------------------------------------------------------|
 * | range               | seqan3::forward_range_concept         | seqan3::view_concept && seqan3::random_access_range_concept |
 * | `range_reference_t` | seqan3::nucleotide_concept            |                                                             |
 *
 * * `irng_t` is the type of the range modified by this view (input).
 * * `rrng_type` is the type of the range returned by this view.
 * * for more details, see \ref view.
 *
 * \par Example
 *
 * Operating on a range of seqan3::dna5:
 * ```cpp
 * dna5_vector vec{"ACGTACGTACGTA"_dna5};
 *
 * // single frame translation
 * auto v1 = vec | view::translate_frames(translation_frames::SINGLE_FRAME);                                              // == [[T,Y,V,R]]
 *
 * // reverse translation
 * auto v2 = vec | view::translate_frames(translation_frames::WITH_REVERSE_COMPLEMENT);                         // == [[T,Y,V,R],[Y,V,R,T]]
 *
 * // forward frames translation
 * auto v3 = vec | view::translate_frames(translation_frames::WITH_FRAME_SHIFTS);                       // == [[T,Y,V,R],[R,T,Y,V],[V,R,T]]
 *
 * // six frame translation
 * auto v4 = vec | view::translate_frames(translation_frames::SIX_FRAME);   // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]

 * // function syntax
 * auto v5 = view::translate_frames(vec, translation_frames::WITH_REVERSE_COMPLEMENT);                          // == [[T,Y,V,R],[Y,V,R,T]]
 *
 * // combinability
 * auto v6 = vec | view::complement | view::translate_frames(translation_frames::WITH_REVERSE_COMPLEMENT);      // == [[C,M,H,A],[M,H,A,C]]
 * ```
 */
seqan3::detail::translate_frames_fn const translate_frames;

} // namespace seqan3::view
