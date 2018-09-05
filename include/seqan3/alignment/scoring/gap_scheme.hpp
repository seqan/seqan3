// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forge, with or without
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
 * \brief Provides seqan3::gap_scheme.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::gap_score
// ------------------------------------------------------------------

/*!\brief A strong type of underlying type `score_type` that represents the score of any character against a gap
 *        character.
 * \tparam score_type The underlying type.
 * \ingroup scoring
 * \see seqan3::gap_scheme
 */
template <arithmetic_concept score_type>
struct gap_score : detail::strong_type<score_type, gap_score<score_type>, detail::strong_type_skill::convert>
{
     using detail::strong_type<score_type, gap_score<score_type>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::gap_score
 * \{
 */
template <arithmetic_concept score_type>
gap_score(score_type &&) -> gap_score<score_type>;
//!\}

// ------------------------------------------------------------------
// seqan3::gap_open_score
// ------------------------------------------------------------------

/*!\brief A strong type of underlying type `score_type` that represents an additional score (usually negative) that
 *        is incurred once additionaly per stretch of consecutive gaps.
 * \tparam score_type The underlying type.
 * \ingroup scoring
 * \see seqan3::gap_scheme
 */
template <arithmetic_concept score_type>
struct gap_open_score : detail::strong_type<score_type, gap_open_score<score_type>, detail::strong_type_skill::convert>
{
     using detail::strong_type<score_type, gap_open_score<score_type>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::gap_open_score
 * \{
 */
template <arithmetic_concept score_type>
gap_open_score(score_type &&) -> gap_open_score<score_type>;
//!\}

// ------------------------------------------------------------------
// seqan3::gap_scheme
// ------------------------------------------------------------------

/*!\brief A scheme for representing and computing scores against gap characters.
 * \tparam score_type Type of the score values saved internally.
 * \ingroup scoring
 */
template <arithmetic_concept score_t = int8_t>
class gap_scheme
{
public:

    /*!\name Member types
     * \{
     */
    //!\brief The template parameter exposed as member type.
    using score_type = score_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr gap_scheme() noexcept = default;
    constexpr gap_scheme(gap_scheme const &) noexcept = default;
    constexpr gap_scheme(gap_scheme &&) noexcept = default;
    constexpr gap_scheme & operator=(gap_scheme const &) noexcept = default;
    constexpr gap_scheme & operator=(gap_scheme &&) noexcept = default;
    ~gap_scheme() noexcept = default;

    /*!\brief Constructor for the Affine gap costs model (delegates to set_affine()).
     * \copydetails set_affine()
     */
    template <arithmetic_concept score_arg_t>
    constexpr gap_scheme(gap_score<score_arg_t> const g, gap_open_score<score_arg_t> const go)
    {
        set_affine(g, go);
    }

    /*!\brief Constructor for the Linear gap costs model (delegates to set_linear()).
     * \copydetails set_linear()
     */
    template <arithmetic_concept score_arg_t>
    constexpr gap_scheme(gap_score<score_arg_t> const g)
    {
        set_linear(g);
    }
    //!\}

    /*!\name Scheme selection
     * \{
     */
    /*!\brief Set the Affine gap costs model.
     * \tparam score_arg_t The underlying type of the arguments.
     * \param[in] g  The cost of each gap character (of type seqan3::gap_score).
     * \param[in] go The additional cost per sequence of gaps (of type seqan3::gap_open_score).
     * \throws std::invalid_argument Thrown if you pass a value that is to large/low to be represented by `score_t`.
     *
     * \details
     *
     * The score for a sequence of `n` gap characters is computed as `n * g + go`.
     *
     * \attention This is the formula used most commonly in the literature, but it is different from SeqAn2 where the
     * formula was `(n-1) * g + go`.
     *
     */
    template <arithmetic_concept score_arg_t>
    constexpr void set_affine(gap_score<score_arg_t> const g, gap_open_score<score_arg_t> const go)
    {
        std::conditional_t<std::Integral<score_t>, int64_t, double> i_g = static_cast<score_arg_t>(g);
        std::conditional_t<std::Integral<score_t>, int64_t, double> i_go = static_cast<score_arg_t>(go);
        if ((i_g  < std::numeric_limits<score_t>::lowest() || i_g  > std::numeric_limits<score_t>::max()) ||
            (i_go < std::numeric_limits<score_t>::lowest() || i_go > std::numeric_limits<score_t>::max()))
        {
            throw std::invalid_argument{"You passed a score value to set_affine/set_linear that is out of range of the "
                                        "scoring scheme's underlying type. Define your scoring scheme with a larger "
                                        "template parameter or down-cast your score value beforehand to prevent "
                                        "this exception."};
        }

        gap = static_cast<score_arg_t>(g);
        gap_open = static_cast<score_arg_t>(go);
    }

    /*!\brief Set the Linear gap costs model.
     * \tparam score_arg_t The underlying type of the argument.
     * \param[in] g The cost of each gap character (of type seqan3::gap_score).
     * \throws std::invalid_argument Thrown if you pass a value that is to large/low to be represented by `score_t`.
     *
     * \details
     *
     * The score for a sequence of `n` gap characters is computed as `n * g`. This is the same as the affine model
     * with a gap open score of `0`.
     */
    template <arithmetic_concept score_arg_t>
    constexpr void set_linear(gap_score<score_arg_t> const g)
    {
        set_affine(g, gap_open_score<score_arg_t>{0});
    }
    //!\}

    /*!\name Accessors
     * \{
     */
    /*!\brief Return the gap score.
     */
    constexpr score_t & get_gap_score() noexcept
    {
        return gap;
    }

    //!\copydoc seqan3::gap_score
    constexpr score_t get_gap_score() const noexcept
    {
        return gap;
    }

    /*!\brief Return the gap open score.
     */
    constexpr score_t & get_gap_open_score() noexcept
    {
        return gap_open;
    }

    //!\copydoc seqan3::gap_open_score
    constexpr score_t get_gap_open_score() const noexcept
    {
        return gap_open;
    }

    /*!\brief Compute the score of a stretch of gap characters.
     * \param number_of_consecutive_gaps The number of consecutive gaps that you wish to know the score for.
     * \returns A signed integer (usually 64bit) that holds the score computed based on the selected scheme.
     */
    constexpr ptrdiff_t score(size_t const number_of_consecutive_gaps) noexcept
    {
        return (gap_open * (number_of_consecutive_gaps ? 1 : 0)) + number_of_consecutive_gaps * gap;
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(gap_scheme const & rhs) const noexcept
    {
        return std::tie(gap, gap_open) == std::tie(rhs.gap, rhs.gap_open);
    }

    constexpr bool operator!=(gap_scheme const & rhs) const noexcept
    {
        return !(*this == rhs);
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive_concept.
     * \param  archive   The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive_concept archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(gap);
        archive(gap_open);
    }
    //!\endcond

private:
    //!\brief The score per gap character.
    score_t gap = 0;
    //!\brief The score per sequence of gaps.
    score_t gap_open = 0;
};

/*!\name Type deduction guides
 * \relates seqan3::gap_scheme
 * \{
 */
gap_scheme() -> gap_scheme<int8_t>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `float`
 * for floating point types.
 * To use a larger type, specify the template argument manually.
 */
template <floating_point_concept score_arg_type>
gap_scheme(gap_score<score_arg_type>, gap_open_score<score_arg_type>) -> gap_scheme<float>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `float`
 * for floating point types.
 * To use a larger type, specify the template argument manually.
 */
template <floating_point_concept score_arg_type>
gap_scheme(gap_score<score_arg_type>) -> gap_scheme<float>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `int8_t`
 * for integer types.
 * To use a larger type, specify the template argument manually.
 */
template <arithmetic_concept score_arg_type>
gap_scheme(gap_score<score_arg_type>, gap_open_score<score_arg_type>) -> gap_scheme<int8_t>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `int8_t`
 * for integer types.
 * To use a larger type, specify the template argument manually.
 */
template <arithmetic_concept score_arg_type>
gap_scheme(gap_score<score_arg_type>) -> gap_scheme<int8_t>;
//!\}

} // namespace seqan3
