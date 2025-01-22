// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides configuration for alignment output.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/utility/concept.hpp>

namespace seqan3::align_cfg
{

/*!\brief Configures the alignment result to output the score.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to compute and output the score. If this option is not set in the alignment
 * configuration, accessing the score via the seqan3::alignment_result object is forbidden and will lead to a compile
 * time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_score.cpp
 *
 * \see seqan3::align_cfg::output_end_position
 * \see seqan3::align_cfg::output_begin_position
 * \see seqan3::align_cfg::output_alignment
 * \see seqan3::align_cfg::output_sequence1_id
 * \see seqan3::align_cfg::output_sequence2_id
 */
class output_score : private pipeable_config_element
{
public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr output_score() = default;                                 //!< Defaulted.
    constexpr output_score(output_score const &) = default;             //!< Defaulted.
    constexpr output_score(output_score &&) = default;                  //!< Defaulted.
    constexpr output_score & operator=(output_score const &) = default; //!< Defaulted.
    constexpr output_score & operator=(output_score &&) = default;      //!< Defaulted.
    ~output_score() = default;                                          //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_score};
};

/*!\brief Configures the alignment result to output the end position.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to compute and output the end positions of the aligned sequences.
 * The end positions must not be identical to the end of the original source sequences. For example,
 * the optimal local alignment might only represent a slice of the original sequences.
 * The end positions denote the end of the alignment within the original sequences, i.e. the positions behind
 * the last aligned characters.
 *
 * If this option is not set in the alignment configuration, then accessing the end positions via the
 * seqan3::alignment_result object is forbidden and will lead to a compile time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_end_position.cpp
 *
 * \see seqan3::align_cfg::output_score
 * \see seqan3::align_cfg::output_begin_position
 * \see seqan3::align_cfg::output_alignment
 * \see seqan3::align_cfg::output_sequence1_id
 * \see seqan3::align_cfg::output_sequence2_id
 */
class output_end_position : private pipeable_config_element
{
public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr output_end_position() = default;                                        //!< Defaulted.
    constexpr output_end_position(output_end_position const &) = default;             //!< Defaulted.
    constexpr output_end_position(output_end_position &&) = default;                  //!< Defaulted.
    constexpr output_end_position & operator=(output_end_position const &) = default; //!< Defaulted.
    constexpr output_end_position & operator=(output_end_position &&) = default;      //!< Defaulted.
    ~output_end_position() = default;                                                 //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_end_position};
};

/*!\brief Configures the alignment result to output the begin positions.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to compute and output the begin positions of the aligned sequences.
 * The begin positions must not be identical to the begin position of the original source sequences. For example,
 * the optimal local alignment might only represent a slice of the original sequences.
 * The begin positions denote the begin of the alignment within the original sequences, i.e. the positions of the
 * first aligned characters.
 *
 * If this option is not set in the alignment configuration, accessing the begin positions via the
 * seqan3::alignment_result object is forbidden and will lead to a compile time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_begin_position.cpp
 *
 * \see seqan3::align_cfg::output_score
 * \see seqan3::align_cfg::output_end_position
 * \see seqan3::align_cfg::output_alignment
 * \see seqan3::align_cfg::output_sequence1_id
 * \see seqan3::align_cfg::output_sequence2_id
 */
class output_begin_position : private pipeable_config_element
{
public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr output_begin_position() = default;                                          //!< Defaulted.
    constexpr output_begin_position(output_begin_position const &) = default;             //!< Defaulted.
    constexpr output_begin_position(output_begin_position &&) = default;                  //!< Defaulted.
    constexpr output_begin_position & operator=(output_begin_position const &) = default; //!< Defaulted.
    constexpr output_begin_position & operator=(output_begin_position &&) = default;      //!< Defaulted.
    ~output_begin_position() = default;                                                   //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_begin_position};
};

/*!\brief Configures the alignment result to output the alignment.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to compute and output the actual aligned sequences.
 *
 * If this option is not set in the alignment configuration, accessing the alignment via the
 * seqan3::alignment_result object is forbidden and will lead to a compile time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_alignment.cpp
 *
 * \see seqan3::align_cfg::output_score
 * \see seqan3::align_cfg::output_end_position
 * \see seqan3::align_cfg::output_begin_position
 * \see seqan3::align_cfg::output_sequence1_id
 * \see seqan3::align_cfg::output_sequence2_id
 */
class output_alignment : private pipeable_config_element
{
public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr output_alignment() = default;                                     //!< Defaulted.
    constexpr output_alignment(output_alignment const &) = default;             //!< Defaulted.
    constexpr output_alignment(output_alignment &&) = default;                  //!< Defaulted.
    constexpr output_alignment & operator=(output_alignment const &) = default; //!< Defaulted.
    constexpr output_alignment & operator=(output_alignment &&) = default;      //!< Defaulted.
    ~output_alignment() = default;                                              //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_alignment};
};

/*!\brief Configures the alignment result to output the id of the first sequence.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to output the id of the first sequence.
 *
 * If this option is not set in the alignment configuration, accessing the id of the first sequence via the
 * seqan3::alignment_result object is forbidden and will lead to a compile time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_sequence1_id.cpp
 *
 * \see seqan3::align_cfg::output_score
 * \see seqan3::align_cfg::output_end_position
 * \see seqan3::align_cfg::output_begin_position
 * \see seqan3::align_cfg::output_alignment
 * \see seqan3::align_cfg::output_sequence2_id
 */
class output_sequence1_id : private pipeable_config_element
{
public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr output_sequence1_id() = default;                                        //!< Defaulted.
    constexpr output_sequence1_id(output_sequence1_id const &) = default;             //!< Defaulted.
    constexpr output_sequence1_id(output_sequence1_id &&) = default;                  //!< Defaulted.
    constexpr output_sequence1_id & operator=(output_sequence1_id const &) = default; //!< Defaulted.
    constexpr output_sequence1_id & operator=(output_sequence1_id &&) = default;      //!< Defaulted.
    ~output_sequence1_id() = default;                                                 //!< Defaulted.

    //!\}
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_sequence1_id};
};

/*!\brief Configures the alignment result to output the id of the second sequence.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This option forces the alignment to output the id of the second sequence.
 *
 * If this option is not set in the alignment configuration, accessing the id of the second sequence via the
 * seqan3::alignment_result object is forbidden and will lead to a compile time error.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_sequence2_id.cpp
 *
 * \see seqan3::align_cfg::output_score
 * \see seqan3::align_cfg::output_end_position
 * \see seqan3::align_cfg::output_begin_position
 * \see seqan3::align_cfg::output_alignment
 * \see seqan3::align_cfg::output_sequence1_id
 */
class output_sequence2_id : private pipeable_config_element
{
public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr output_sequence2_id() = default;                                        //!< Defaulted.
    constexpr output_sequence2_id(output_sequence2_id const &) = default;             //!< Defaulted.
    constexpr output_sequence2_id(output_sequence2_id &&) = default;                  //!< Defaulted.
    constexpr output_sequence2_id & operator=(output_sequence2_id const &) = default; //!< Defaulted.
    constexpr output_sequence2_id & operator=(output_sequence2_id &&) = default;      //!< Defaulted.
    ~output_sequence2_id() = default;                                                 //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::output_sequence2_id};
};

} // namespace seqan3::align_cfg
