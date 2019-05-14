// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::simd_gap_scheme.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/scoring/gap_scheme_base.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::simd_gap_scheme
// ------------------------------------------------------------------

/*!\brief A scheme for representing and computing scores against gap characters.
 * \tparam score_type Type of the score values saved internally.
 * \ingroup scoring
 */
template <simd_concept simd_t>
class simd_gap_scheme : public gap_scheme_base<simd_gap_scheme<simd_t>, simd_t>
{
private:

    //!\brief The type of the base class.
    using base_t = gap_scheme_base<simd_gap_scheme<simd_t>, simd_t>;

    //!\befriend the base_t.
    friend base_t;

public:

    using base_t::score_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr simd_gap_scheme() noexcept = default;                                    //!< Defaulted
    constexpr simd_gap_scheme(simd_gap_scheme const &) noexcept = default;             //!< Defaulted
    constexpr simd_gap_scheme(simd_gap_scheme &&) noexcept = default;                  //!< Defaulted
    constexpr simd_gap_scheme & operator=(simd_gap_scheme const &) noexcept = default; //!< Defaulted
    constexpr simd_gap_scheme & operator=(simd_gap_scheme &&) noexcept = default;      //!< Defaulted
    ~simd_gap_scheme() noexcept = default;                                             //!< Defaulted

    /*!\brief Constructor for the Affine gap costs model (delegates to set_affine()).
     */
    template <Arithmetic score_arg_t>
    //!\cond
        requires std::ConvertibleTo<typename simd_traits<simd_t>::scalar_type, score_arg_t>
    //!\endcond
    constexpr simd_gap_scheme(gap_score<score_arg_t> const g, gap_open_score<score_arg_t> const go) : base_t{g, go}
    {}

    /*!\brief Constructor for the Linear gap costs model (delegates to set_linear()).
     */
    template <Arithmetic score_arg_t>
    //!\cond
        requires std::ConvertibleTo<typename simd_traits<simd_t>::scalar_type, score_arg_t>
    //!\endcond
    constexpr simd_gap_scheme(gap_score<score_arg_t> const g) : base_t{g}
    {}
    //!\}

private:

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
    template <Arithmetic score_arg_t>
    //!\cond
        requires std::ConvertibleTo<typename simd_traits<simd_t>::scalar_type, score_arg_t>
    //!\endcond
    constexpr void set_scheme_impl(score_arg_t const g, score_arg_t const go)
    {
        using score_t = typename simd_traits<simd_t>::scalar_type;
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

        base_t::gap = simd::fill<simd_t>(g);
        base_t::gap_open = simd::fill<simd_t>(go);
    }
    //!\}

    /*!\name Accessors
     * \{
     */

    /*!\brief Compute the score of a stretch of gap characters.
     * \param number_of_consecutive_gaps The number of consecutive gaps that you wish to know the score for.
     * \returns A signed integer (usually 64bit) that holds the score computed based on the selected scheme.
     */
    constexpr ptrdiff_t score_impl(size_t const number_of_consecutive_gaps) const noexcept
    {
        return (base_t::gap_open[0] * (number_of_consecutive_gaps ? 1 : 0)) +
                number_of_consecutive_gaps * base_t::gap[0];
    }
    //!\}
};

} // namespace seqan3
