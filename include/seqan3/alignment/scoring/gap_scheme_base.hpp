// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::gap_scheme_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
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
template <Arithmetic score_type>
struct gap_score : detail::strong_type<score_type, gap_score<score_type>, detail::strong_type_skill::convert>
{
     using detail::strong_type<score_type, gap_score<score_type>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::gap_score
 * \{
 */
template <Arithmetic score_type>
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
template <Arithmetic score_type>
struct gap_open_score : detail::strong_type<score_type, gap_open_score<score_type>, detail::strong_type_skill::convert>
{
     using detail::strong_type<score_type, gap_open_score<score_type>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::gap_open_score
 * \{
 */
template <Arithmetic score_type>
gap_open_score(score_type &&) -> gap_open_score<score_type>;
//!\}

template <typename derived_t, typename score_t>
class gap_scheme_base
{
public:

    /*!\name Member types
     * \{
     */
    //!\brief The template parameter exposed as member type.
    using score_type = score_t;
    //!\}

private:

    //!\brief Befriend the derived type.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr gap_scheme_base(gap_scheme_base const &) noexcept = default;             //!< Default.
    constexpr gap_scheme_base(gap_scheme_base &&) noexcept = default;                  //!< Default.
    constexpr gap_scheme_base & operator=(gap_scheme_base const &) noexcept = default; //!< Default.
    constexpr gap_scheme_base & operator=(gap_scheme_base &&) noexcept = default;      //!< Default.
    ~gap_scheme_base() noexcept = default;                                             //!< Default.

    constexpr gap_scheme_base() noexcept
    {
        derived().set_scheme_impl(-1, 0);
    }

    /*!\brief Constructor for the Affine gap costs model (delegates to set_affine()).
     * \copydetails set_affine()
     */
    template <Arithmetic score_arg_t>
    constexpr gap_scheme_base(gap_score<score_arg_t> const g, gap_open_score<score_arg_t> const go)
    {
        set_affine(g, go);
    }

    /*!\brief Constructor for the Linear gap costs model (delegates to set_linear()).
     * \copydetails set_linear()
     */
    template <Arithmetic score_arg_t>
    constexpr gap_scheme_base(gap_score<score_arg_t> const g)
    {
        set_linear(g);
    }
    //!\}

public:

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
    constexpr void set_affine(gap_score<score_arg_t> const g, gap_open_score<score_arg_t> const go)
    {
        derived().set_scheme_impl(g.get(), go.get());
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
    template <Arithmetic score_arg_t>
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
    constexpr score_t const & get_gap_score() const noexcept
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
    constexpr score_t const & get_gap_open_score() const noexcept
    {
        return gap_open;
    }

    /*!\brief Compute the score of a stretch of gap characters.
     * \param number_of_consecutive_gaps The number of consecutive gaps that you wish to know the score for.
     * \returns A signed integer (usually 64bit) that holds the score computed based on the selected scheme.
     */
    constexpr ptrdiff_t score(size_t const number_of_consecutive_gaps) const noexcept
    {
        return derived().score_impl(number_of_consecutive_gaps);
        // return (gap_open * (number_of_consecutive_gaps ? 1 : 0)) + number_of_consecutive_gaps * gap;
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(gap_scheme_base const & rhs) const noexcept
    {
        return std::tie(gap, gap_open) ==
               std::tie(rhs.gap, rhs.gap_open);
    }

    constexpr bool operator!=(gap_scheme_base const & rhs) const noexcept
    {
        return !(*this == rhs);
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::CerealArchive.
     * \param  archive   The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <CerealArchive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(get_gap_score());
        archive(get_gap_open_score());
    }
    //!\endcond

private:

    derived_t & derived() noexcept
    {
        return static_cast<derived_t &>(*this);
    }

    derived_t const & derived() const noexcept
    {
        return static_cast<derived_t const &>(*this);
    }

    //!\brief The score per gap character. Defaults to -1.
    score_t gap{};
    //!\brief The score per sequence of gaps. Defaults to 0.
    score_t gap_open{};
};

} // namespace seqan3
