// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the configuration for maximum number of errors for all error types.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <variant>

#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/configuration/max_error_common.hpp>
#include <seqan3/search/detail/search_common.hpp>

namespace seqan3::search_cfg
{
/*!\brief Configuration element that represents the number or rate of total errors.
 * \ingroup search_configuration
 * \see search_configuration
 *
 * \details
 * This configuration element can be used to determine the number or rate of total errors that are supported.
 *
 * ### Example
 * \include test/snippet/search/configuration_error.cpp
 */
class max_error_total : private pipeable_config_element
{
public:
    //!\brief The error count or error rate.
    std::variant<error_count, error_rate> error{error_count{0}};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    max_error_total() = default;                                    //!< Defaulted.
    max_error_total(max_error_total const &) = default;             //!< Defaulted.
    max_error_total(max_error_total &&) = default;                  //!< Defaulted.
    max_error_total & operator=(max_error_total const &) = default; //!< Defaulted.
    max_error_total & operator=(max_error_total &&) = default;      //!< Defaulted.
    ~max_error_total() = default;                                   //!< Defaulted.

    /*!\brief Initialises the total error with the given seqan3::search_cfg::error_count.
     * \param[in] error The maximal number of total errors allowed in the search.
     */
    constexpr explicit max_error_total(error_count error) : error{std::move(error)}
    {}

    /*!\brief Initialises the total error with the given seqan3::search_cfg::error_rate.
     * \param[in] error The maximal total error rate allowed in the search.
     */
    constexpr explicit max_error_total(error_rate error) : error{std::move(error)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::max_error_total};
};

/*!\brief Configuration element that represents the number or rate of substitution errors.
 * \ingroup search_configuration
 * \see search_configuration
 *
 * \details
 * This configuration element can be used to determine the number or rate of substitution errors that are supported.
 * A substitution corresponds to diverging bases between text and query for a certain position.
 *
 * ### Example
 * \include test/snippet/search/configuration_error.cpp
 */
class max_error_substitution : private pipeable_config_element
{
public:
    //!\brief The error count or error rate.
    std::variant<error_count, error_rate> error{error_count{0}};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    max_error_substitution() = default;                                           //!< Defaulted.
    max_error_substitution(max_error_substitution const &) = default;             //!< Defaulted.
    max_error_substitution(max_error_substitution &&) = default;                  //!< Defaulted.
    max_error_substitution & operator=(max_error_substitution const &) = default; //!< Defaulted.
    max_error_substitution & operator=(max_error_substitution &&) = default;      //!< Defaulted.
    ~max_error_substitution() = default;                                          //!< Defaulted.

    /*!\brief Initialises the substitution error with the given seqan3::search_cfg::error_count.
     * \param[in] error The maximal number of substitution errors allowed in the search.
     */
    constexpr explicit max_error_substitution(error_count error) : error{std::move(error)}
    {}

    /*!\brief Initialises the substitution error with the given seqan3::search_cfg::error_rate.
     * \param[in] error The maximal error rate for substitutions allowed in the search.
     */
    constexpr explicit max_error_substitution(error_rate error) : error{std::move(error)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::max_error_substitution};
};

/*!\brief Configuration element that represents the number or rate of insertion errors.
 * \ingroup search_configuration
 * \see search_configuration
 *
 * \details
 * This configuration element can be used to determine the number or rate of insertion errors that are supported.
 * An insertion corresponds to a base inserted into the query that does not occur in the text at the position.
 *
 * ### Example
 * \include test/snippet/search/configuration_error.cpp
 */
class max_error_insertion : private pipeable_config_element
{
public:
    //!\brief The error count or error rate.
    std::variant<error_count, error_rate> error{error_count{0}};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    max_error_insertion() = default;                                        //!< Defaulted.
    max_error_insertion(max_error_insertion const &) = default;             //!< Defaulted.
    max_error_insertion(max_error_insertion &&) = default;                  //!< Defaulted.
    max_error_insertion & operator=(max_error_insertion const &) = default; //!< Defaulted.
    max_error_insertion & operator=(max_error_insertion &&) = default;      //!< Defaulted.
    ~max_error_insertion() = default;                                       //!< Defaulted.

    /*!\brief Initialises the insertion error with the given seqan3::search_cfg::error_count.
     * \param[in] error The maximal number of insertion errors allowed in the search.
     */
    constexpr explicit max_error_insertion(error_count error) : error{std::move(error)}
    {}

    /*!\brief Initialises the insertion error with the given seqan3::search_cfg::error_rate.
     * \param[in] error The maximal error rate for insertions allowed in the search.
     */
    constexpr explicit max_error_insertion(error_rate error) : error{std::move(error)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::max_error_insertion};
};

/*!\brief Configuration element that represents the number or rate of deletion errors.
 * \ingroup search_configuration
 * \see search_configuration
 *
 * \details
 * This configuration element can be used to determine the number or rate of deletion errors that are supported.
 * A deletion corresponds to a base deleted from the query sequence that does occur in the text.
 * Deletions at the beginning and at the end of the sequence are not considered during a search.
 *
 * ### Example
 * \include test/snippet/search/configuration_error.cpp
 */
class max_error_deletion : private pipeable_config_element
{
public:
    //!\brief The error count or error rate.
    std::variant<error_count, error_rate> error{error_count{0}};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    max_error_deletion() = default;                                       //!< Defaulted.
    max_error_deletion(max_error_deletion const &) = default;             //!< Defaulted.
    max_error_deletion(max_error_deletion &&) = default;                  //!< Defaulted.
    max_error_deletion & operator=(max_error_deletion const &) = default; //!< Defaulted.
    max_error_deletion & operator=(max_error_deletion &&) = default;      //!< Defaulted.
    ~max_error_deletion() = default;                                      //!< Defaulted.

    /*!\brief Initialises the deletion error with the given seqan3::search_cfg::error_count.
     * \param[in] error The maximal number of deletion errors allowed in the search.
     */
    constexpr explicit max_error_deletion(error_count error) : error{std::move(error)}
    {}

    /*!\brief Initialises the deletion error with the given seqan3::search_cfg::error_rate.
     * \param[in] error The maximal error rate for deletions allowed in the search.
     */
    constexpr explicit max_error_deletion(error_rate error) : error{std::move(error)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::max_error_deletion};
};

} // namespace seqan3::search_cfg
