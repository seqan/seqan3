// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configuration for alignment output.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/concept/core_language.hpp>

namespace seqan3::detail
{

/*!\brief Triggers score-only computation of the sequence alignment.
 * \ingroup alignment_configuration
 */
struct with_score_type
{
    //!\privatesection
    //!\brief An internal rank used for an ordered access of seqan3::align_cfg::result options.
    static constexpr int8_t rank = 0;
};

/*!\brief Triggers score computation and determines the end position of the sequence alignment.
 * \ingroup alignment_configuration
 */
struct with_back_coordinate_type
{
    //!\privatesection
    //!\brief An internal rank used for an ordered access of seqan3::align_cfg::result options.
    static constexpr int8_t rank = 1;
};

/*!\brief Triggers score computation and determines begin and end position of the sequence alignment.
 * \ingroup alignment_configuration
 */
struct with_front_coordinate_type
{
    //!\privatesection
    //!\brief An internal rank used for an ordered access of seqan3::align_cfg::result options.
    static constexpr int8_t rank = 2;
};

/*!\brief Triggers score computation and determines the end position of the sequence alignment as well as the
 *        full traceback.
 * \ingroup alignment_configuration
 */
struct with_alignment_type
{
    //!\privatesection
    //!\brief An internal rank used for an ordered access of seqan3::align_cfg::result options.
    static constexpr int8_t rank = 3;
};

/*!\brief Helper type to configure the score type of the alignment algorithm.
 * \ingroup alignment_configuration
 */
template <typename t>
struct score_type
{};

} // namespace seqan3::detail

namespace seqan3
{

//!\brief Helper variable used to select score-only computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_score_type with_score{};
//!\brief Helper variable used to select end-position computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_back_coordinate_type with_back_coordinate{};
//!\brief Helper variable used to select begin position computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_front_coordinate_type with_front_coordinate{};
//!\brief Helper Variable used to select trace computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_alignment_type with_alignment{};
/*!\brief Helper variable used to configure the score type for the alignment algorithm.
 * \relates seqan3::align_cfg::result
 * \tparam t The type to use for the computed alignment score; must model seqan3::arithmetic.
 */
template <arithmetic score_t>
inline constexpr detail::score_type<score_t> using_score_type{};

} // namespace seqan3

namespace seqan3::align_cfg
{

/*!\brief Sets the result of the alignment computation.
 * \ingroup alignment_configuration
 * \tparam alignment_result_tag_t The type used to specify which feature should be computed during the pairwise
 *                                alignment. Defaults to seqan3::detail::with_score_type.
 *
 * \details
 *
 * The output of the pairwise alignment can be configured using this result configuration element. Depending on the
 * settings, the most efficient implementation is chosen to compute the result.
 * Currently four different modes can be configured (first constructor parameter):
 *
 * 1. computing only the \ref seqan3::align_cfg::result::with_score "score",
 * 2. computing in addition the \ref seqan3::align_cfg::result::with_back_coordinate "end position",
 * 3. computing in addition the \ref seqan3::align_cfg::result::with_front_coordinate "begin position",
 * 4. and finally also computing the \ref seqan3::align_cfg::result::with_alignment "alignment".
 *
 * These settings will directly affect the contents of the seqan3::alignment_result object which is returned by the
 * alignment algorithm. For example, if you chose the \ref seqan3::align_cfg::result::with_alignment "alignment"
 * feature, your result object will contain the score, end point, begin point and the alignment.
 *
 * In addition, you can specify the \ref seqan3::align_cfg::result::using_score_type "score type"
 * with the second constructor argument (see example).
 *
 * By default, the alignment algorithm will only compute the score with score type `int32_t`.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_result_example.cpp
 */
template <typename alignment_result_tag_t = detail::with_score_type, typename score_t = int32_t>
//!\cond
    requires std::same_as<alignment_result_tag_t, detail::with_score_type> ||
             std::same_as<alignment_result_tag_t, detail::with_back_coordinate_type> ||
             std::same_as<alignment_result_tag_t, detail::with_front_coordinate_type> ||
             std::same_as<alignment_result_tag_t, detail::with_alignment_type>
//!\endcond
class result : public pipeable_config_element<result<alignment_result_tag_t, score_t>, alignment_result_tag_t>
{
    //!\brief The base type of this class.
    using base_t = pipeable_config_element<result<alignment_result_tag_t, score_t>, alignment_result_tag_t>;
public:
    //!\brief The score type of the alignment result.
    using score_type = score_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr result() = default; //!< Defaulted.
    constexpr result(result const &) = default; //!< Defaulted.
    constexpr result & operator=(result const &) = default; //!< Defaulted.
    constexpr result(result &&) = default; //!< Defaulted.
    constexpr result & operator=(result &&) = default; //!< Defaulted.
    ~result() = default; //!< Defaulted.

    /*!\brief Construction from the result feature you want to compute (e.g. seqan3::with_score).
     * \param[in] result_tag The feature you want the alignment algorithm to compute (e.g. seqan3::with_score).
     */
    constexpr result(alignment_result_tag_t result_tag) noexcept : base_t{result_tag} {}

    /*!\brief Construction from the result feature you want to compute (e.g. seqan3::with_score).
     * \param[in] result_tag The feature you want the alignment algorithm to compute (e.g. seqan3::with_score).
     * \param[in] score_type_tag The score type to use (e.g. seqan3::using_score_type<int>).
     */
    constexpr result(alignment_result_tag_t result_tag,
                     detail::score_type<score_t> SEQAN3_DOXYGEN_ONLY(score_type_tag)) noexcept : base_t{result_tag} {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::result};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::result
 * \{
 */
//!\brief Deduces the alignment result from the given constructor argument.
template <typename alignment_result_tag_t>
result(alignment_result_tag_t) -> result<alignment_result_tag_t>;

//!\brief Deduces the alignment result from the given constructor arguments.
template <typename alignment_result_tag_t, arithmetic score_t>
result(alignment_result_tag_t, detail::score_type<score_t>) -> result<alignment_result_tag_t, score_t>;
//!\}
} // namespace seqan3::align_cfg
