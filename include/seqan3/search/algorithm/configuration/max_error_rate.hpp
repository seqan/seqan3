// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
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
 * \brief Provides the configuration for maximum number of errors in percent of the query length across all error types.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/algorithm/fill.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/slice.hpp>

#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/search/algorithm/configuration/detail.hpp>
#include <seqan3/search/algorithm/configuration/max_error_common.hpp>

namespace seqan3::search_cfg
{

/*!\brief A configuration element for the maximum number of errors in percent of the query length across all error types
 *        (mismatches, insertions, deletions). This is an upper bound of errors independent from error rates of
 *        specific error types.
 * \details An insertion corresponds to a base inserted into the query that does not occur in the text at the position,
 *          a deletion corresponds to a base deleted from the query sequence that does occur in the indexed text.
 *          Deletions at the beginning and at the end of the sequence are not considered during a search.
 * \ingroup search_configuration
 */
template <typename ...errors_t>
//!\cond
    requires sizeof...(errors_t) <= 4 &&
            ((detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, total> ||
              detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, substitution> ||
              detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, deletion>  ||
              detail::is_type_specialisation_of_v<std::remove_reference_t<errors_t>, insertion>) && ...)
//!\endcond
class max_error_rate : public pipeable_config_element
{

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

    //!\publicsecton
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr max_error_rate()                                   noexcept = default;
    constexpr max_error_rate(max_error_rate const &)             noexcept = default;
    constexpr max_error_rate(max_error_rate &&)                  noexcept = default;
    constexpr max_error_rate & operator=(max_error_rate const &) noexcept = default;
    constexpr max_error_rate & operator=(max_error_rate &&)      noexcept = default;
    ~max_error_rate()                                            noexcept = default;

    /*!\brief Constructs the object from a set of error specifiers.
     * \tparam    errors_t A template parameter pack with the error types.
     * \param[in] errors   A pack of error specifiers.
     *
     * \details
     *
     * \todo write me
     */
    constexpr max_error_rate(errors_t && ...errors)
    //!\cond
        requires sizeof...(errors_t) > 0
    //!\endcond
    {
        detail::for_each_value([this](auto e)
        {
            value[remove_cvref_t<decltype(e)>::_id()] = e.get();
        }, std::forward<errors_t>(errors)...);

        // check correct values.
        ranges::for_each(value, [](auto error_elem)
        {
            if (0.0 > error_elem  || error_elem > 1.0)
                throw std::invalid_argument("Error rates must be between 0 and 1.");
        });

        // Only total is set so we set all other errors to the total limit.
        if constexpr (((std::remove_reference_t<errors_t>::_id() == 0) || ...) && sizeof...(errors) == 1)
        {
            ranges::fill(value | ranges::view::slice(1, 4), value[0]);
        } // otherwise if total is not set but any other field is set than use total as the sum of all set errors.
        else if constexpr (!((std::remove_reference_t<errors_t>::_id() == 0) || ...) && sizeof...(errors) > 0)
        {
            value[0] = std::min(1., ranges::accumulate(value | ranges::view::slice(1, 4), .0));
        }
    }
    //!}

    //!\brief The ordered error values.
    std::array<double, 4> value{.0, .0, .0, .0};
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
