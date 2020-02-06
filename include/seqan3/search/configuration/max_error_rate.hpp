// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the configuration for maximum number of errors in percent of the query length across all error types.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/numeric/accumulate.hpp>

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/configuration/max_error_common.hpp>
#include <seqan3/std/algorithm>

namespace seqan3::search_cfg
{

/*!\brief A configuration element for the maximum number of errors in percent of the query length across all error types
 *        (mismatches, insertions, deletions). This is an upper bound of errors independent from error rates of
 *        specific error types.
 * \ingroup search_configuration
 * \details A mismatch corresponds to diverging bases between text and query for a certain position.
 *          An insertion corresponds to a base inserted into the query that does not occur in the text at the position,
 *          a deletion corresponds to a base deleted from the query sequence that does occur in the indexed text.
 *          Deletions at the beginning and at the end of the sequence are not considered during a search.
 */
template <typename ...errors_t>
//!\cond
    requires sizeof...(errors_t) <= 4 &&
            ((detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, total> ||
              detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, substitution> ||
              detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, deletion>  ||
              detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, insertion>) && ...)
//!\endcond
class max_error_rate : public pipeable_config_element<max_error_rate<errors_t...>, std::array<double, 4>>
{

    //!\brief An alias type for the base class.
    using base_t = pipeable_config_element<max_error_rate<errors_t...>, std::array<double, 4>>;

    //!\brief Helper function to check valid error rate configuration.
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

    static_assert(check_consistency(errors_t{}...),
                  "You may not use the same error specifier more than once.");

public:

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::max_error_rate};

    //!\publicsection
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr max_error_rate()                                   noexcept = default; //!< Default constructor.
    constexpr max_error_rate(max_error_rate const &)             noexcept = default; //!< Copy constructor.
    constexpr max_error_rate(max_error_rate &&)                  noexcept = default; //!< Move constructor.
    constexpr max_error_rate & operator=(max_error_rate const &) noexcept = default; //!< Copy assignment.
    constexpr max_error_rate & operator=(max_error_rate &&)      noexcept = default; //!< Move assignment.
    ~max_error_rate()                                            noexcept = default; //!< Destructor.

    /*!\brief Constructs the object from a set of error specifiers.
     * \tparam    errors_t A template parameter pack with the error types.
     * \param[in] errors   A pack of error specifiers.
     *
     * \details
     *
     * \details
     *
     * This configuration can be used to specify the total rates of error types.
     * It restricts the rates of substitutions, insertions, deletions and total errors within the search to the given
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
     * \include test/snippet/search/configuration_error_rate.cpp
     */
    constexpr max_error_rate(errors_t && ...errors)
    //!\cond
        requires sizeof...(errors_t) > 0
    //!\endcond
        : base_t{}
    {
        detail::for_each([this](auto e)
        {
            base_t::value[remove_cvref_t<decltype(e)>::_id()] = e.get();
        }, std::forward<errors_t>(errors)...);

        // check correct values.
        std::ranges::for_each(base_t::value, [](auto error_elem)
        {
            if (0.0 > error_elem  || error_elem > 1.0)
                throw std::invalid_argument("Error rates must be between 0 and 1.");
        });

        // Only total is set so we set all other errors to the total limit.
        if constexpr (((std::remove_reference_t<errors_t>::_id() == 0) || ...) && sizeof...(errors) == 1)
        {
            std::ranges::fill(base_t::value | views::slice(1, 4), base_t::value[0]);
        } // otherwise if total is not set but any other field is set than use total as the sum of all set errors.
        else if constexpr (!((std::remove_reference_t<errors_t>::_id() == 0) || ...) && sizeof...(errors) > 0)
        {
            base_t::value[0] = std::min(1., ranges::accumulate(base_t::value | views::slice(1, 4), .0));
        }
    }
    //!}
};

/*!\name Type deduction guides
 * \relates seqan3::search_cfg::max_error_rate
 * \{
 */

//!\brief Deduces empty list of error specifiers.
max_error_rate() -> max_error_rate<>;

//!\brief Deduces template arguments from the passed error specifiers.
template <typename ...errors_t>
max_error_rate(errors_t && ...) -> max_error_rate<remove_cvref_t<errors_t>...>;
//!\}

} // namespace seqan3::search_cfg
