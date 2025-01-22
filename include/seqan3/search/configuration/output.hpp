// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the configuration for the content of the search result.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/search/configuration/detail.hpp>

namespace seqan3::search_cfg
{

/*!\brief Include the query_id in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
class output_query_id : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr output_query_id() = default;                                    //!< Defaulted.
    constexpr output_query_id(output_query_id const &) = default;             //!< Defaulted.
    constexpr output_query_id(output_query_id &&) = default;                  //!< Defaulted.
    constexpr output_query_id & operator=(output_query_id const &) = default; //!< Defaulted.
    constexpr output_query_id & operator=(output_query_id &&) = default;      //!< Defaulted.
    ~output_query_id() = default;                                             //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_query_id};
};

/*!\brief Include the reference_id in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
class output_reference_id : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr output_reference_id() = default;                                        //!< Defaulted.
    constexpr output_reference_id(output_reference_id const &) = default;             //!< Defaulted.
    constexpr output_reference_id(output_reference_id &&) = default;                  //!< Defaulted.
    constexpr output_reference_id & operator=(output_reference_id const &) = default; //!< Defaulted.
    constexpr output_reference_id & operator=(output_reference_id &&) = default;      //!< Defaulted.
    ~output_reference_id() = default;                                                 //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_reference_id};
};

/*!\brief Include the reference_begin_position in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
class output_reference_begin_position : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr output_reference_begin_position() = default;                                        //!< Defaulted.
    constexpr output_reference_begin_position(output_reference_begin_position const &) = default; //!< Defaulted.
    constexpr output_reference_begin_position(output_reference_begin_position &&) = default;      //!< Defaulted.
    constexpr output_reference_begin_position &
    operator=(output_reference_begin_position const &) = default;                                        //!< Defaulted.
    constexpr output_reference_begin_position & operator=(output_reference_begin_position &&) = default; //!< Defaulted.
    ~output_reference_begin_position() = default;                                                        //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_reference_begin_position};
};

/*!\brief Include the index_cursor in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
class output_index_cursor : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr output_index_cursor() = default;                                        //!< Defaulted.
    constexpr output_index_cursor(output_index_cursor const &) = default;             //!< Defaulted.
    constexpr output_index_cursor(output_index_cursor &&) = default;                  //!< Defaulted.
    constexpr output_index_cursor & operator=(output_index_cursor const &) = default; //!< Defaulted.
    constexpr output_index_cursor & operator=(output_index_cursor &&) = default;      //!< Defaulted.
    ~output_index_cursor() = default;                                                 //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_index_cursor};
};

} // namespace seqan3::search_cfg
