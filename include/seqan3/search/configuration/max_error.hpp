// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the configuration for maximum number of errors for all error types.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/configuration/max_error_common.hpp>
#include <seqan3/search/detail/search_common.hpp>

namespace seqan3::search_cfg
{
/*!\brief A configuration element for the maximum number of errors across all error types (mismatches, insertions,
 *        deletions). This is an upper bound of errors independent from error numbers of specific error types.
 * \ingroup search_configuration
 * \details A mismatch corresponds to diverging bases between text and query for a certain position.
 *          An insertion corresponds to a base inserted into the query that does not occur in the text at the position,
 *          a deletion corresponds to a base deleted from the query sequence that does occur in the indexed text.
 *          Deletions at the beginning and at the end of the sequence are not considered during a search.
 */
class max_error : public pipeable_config_element<max_error, detail::search_param>
{
    //!\brief An alias type for the base class.
    using base_t = pipeable_config_element<max_error, detail::search_param>;

    //!\brief Helper function to check valid max error configuration.
    template <typename ..._errors_t>
    static constexpr bool check_consistency(_errors_t ...errors)
    {
        if constexpr (sizeof...(errors) < 2)
        {
            return true;
        }
        else
        {
            return [] (auto head, auto ...tail) constexpr
            {
                using head_t = decltype(head);
                if constexpr (((head_t::_id != decltype(tail)::_id) && ...))
                    return check_consistency(tail...);
                else
                    return false;
            }(errors...);
        }
    }

public:
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::max_error};

    //!\publicsection
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr max_error() = default; //!< Defaulted.
    constexpr max_error(max_error const &) = default; //!< Defaulted.
    constexpr max_error(max_error &&) = default; //!< Defaulted.
    constexpr max_error & operator=(max_error const &) = default; //!< Defaulted.
    constexpr max_error & operator=(max_error &&) = default; //!< Defaulted.
    ~max_error() = default; //!< Defaulted

    /*!\brief Constructs the object from a set of error specifiers.
     * \tparam    errors_t A template parameter pack with the error types.
     * \param[in] errors   A pack of error specifiers.
     *
     * \details
     *
     * This configuration can be used to specify the total number of error types.
     * It restricts the number of substitutions, insertions, deletions and total errors within the search to the given
     * values and will behave as follows:
     * |                                  |                              |                           |                        |                       |                                       |                                      |                                   |
     * |----------------------------------|:----------------------------:|:-------------------------:|:----------------------:|:---------------------:|:-------------------------------------:|:------------------------------------:|:---------------------------------:|
     * | **Behaviour**                    | Set all error types to total | Set total to substitution | Set total to insertion | Set total to deletion | Set total to substitution + insertion | Set total to substitution + deletion | Set total to insertion + deletion |
     * | seqan3::search_cfg::total        |               ✅              |                           |                        |                       |                                       |                                      |                                   |
     * | seqan3::search_cfg::substitution |                              |             ✅             |                        |                       |                   ✅                   |                   ✅                  |                                   |
     * | seqan3::search_cfg::insertion    |                              |                           |            ✅           |                       |                   ✅                   |                                      |                 ✅                 |
     * | seqan3::search_cfg::deletion     |                              |                           |                        |           ✅           |                                       |                   ✅                  |                 ✅                 |
     * If seqan3::search_cfg::total and any other error type are specified, all types are set to the respective values.
     *
     * ### Example
     *
     * \include test/snippet/search/configuration_error.cpp
     */
    template <typename ...errors_t>
    //!\cond
        requires sizeof...(errors_t) > 0 && sizeof...(errors_t) <= 4 &&
                 ((detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, total> ||
                   detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, substitution> ||
                   detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, deletion>  ||
                   detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, insertion>) && ...)
    //!\endcond
    constexpr max_error(errors_t && ...errors) noexcept : base_t{}
    {
        static_assert(check_consistency(errors_t{}...), "You may not use the same error specifier more than once.");

        detail::for_each([this](auto e)
        {
            switch (remove_cvref_t<decltype(e)>::_id())
            {
                case 0 : value.total = e.get(); break;
                case 1 : value.substitution = e.get(); break;
                case 2 : value.insertion = e.get(); break;
                case 3 : value.deletion = e.get(); break;
            }
        }, std::forward<errors_t>(errors)...);

        // Only total is set so we set all other errors to the total limit.
        if constexpr (((std::remove_reference_t<errors_t>::_id() == 0) || ...) && sizeof...(errors) == 1)
        {
            value.substitution = value.insertion = value.deletion = value.total;
        } // otherwise if total is not set but any other field is set than use total as the sum of all set errors.
        else if constexpr (!((std::remove_reference_t<errors_t>::_id() == 0) || ...) && sizeof...(errors) > 0)
        {
            value.total = std::min<uint32_t>(255, value.substitution + value.insertion + value.deletion);
        }
    }
    //!\}
};

} // namespace seqan3::search_cfg
