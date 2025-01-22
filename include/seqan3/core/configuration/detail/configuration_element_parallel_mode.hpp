// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::parallel_mode.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <optional>

#include <seqan3/core/configuration/pipeable_config_element.hpp>

namespace seqan3::detail
{
/*!\brief A global configuration type used to enable parallel execution of algorithms.
 * \ingroup core_configuration
 * \tparam wrapped_config_id_t The algorithm specific configuration id wrapped in a std::integral_constant.
 *
 * \details
 *
 * This type is used to enable the parallel mode of the algorithms.
 */
template <typename wrapped_config_id_t>
class parallel_mode : private pipeable_config_element
{
public:
    /*!\name Constructors, assignment and destructor
     * \{
     */
    parallel_mode() = default;                                  //!< Defaulted.
    parallel_mode(parallel_mode const &) = default;             //!< Defaulted.
    parallel_mode(parallel_mode &&) = default;                  //!< Defaulted.
    parallel_mode & operator=(parallel_mode const &) = default; //!< Defaulted.
    parallel_mode & operator=(parallel_mode &&) = default;      //!< Defaulted.
    ~parallel_mode() = default;                                 //!< Defaulted.

    /*!\brief Sets the number of threads for the parallel configuration element.
     * \param[in] thread_count_ The maximum number of threads to be used by the algorithm.
     */
    explicit parallel_mode(uint32_t thread_count_) noexcept : thread_count{thread_count_}
    {}
    //!\}

    //!\brief The maximum number of threads the algorithm can use.
    std::optional<uint32_t> thread_count{std::nullopt};

    /*!\privatesection
     * \brief Internal id to check for consistent configuration settings.
     */
    static constexpr typename wrapped_config_id_t::value_type id{wrapped_config_id_t::value};
};
} // namespace seqan3::detail
