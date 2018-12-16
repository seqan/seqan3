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
 * \brief Provides the configuration for maximum number of errors for all error types.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/search/algorithm/configuration/max_error_common.hpp>
#include <seqan3/search/algorithm/configuration/utility.hpp>

namespace seqan3::detail
{
/*!\brief A configuration element for the maximum number of errors across all error types (mismatches, insertions,
          deletions). This is an upper bound of errors independent from error numbers or rates of specific error types.
 * \ingroup search_configuration
 */
struct search_config_max_error
{
    //!\brief The actual value.
    std::tuple<uint8_t, uint8_t, uint8_t, uint8_t> value;
};

/*!\brief The max_error adaptor enabling pipe notation.
 * \ingroup search_configuration
 */
struct search_config_max_error_adaptor : public configuration_fn_base<search_config_max_error_adaptor>
{
private:
    /*!\brief Helper function to set single error types.
     * \tparam error_type_t Error type to be extracted from `in` and set in `out`.
     * \param[in] out Output tuple containing all error types.
     * \param[in] in Input tuple containing a subset of error types.
     */
    template <typename error_type_t, typename out_t, typename ... error_types_t>
    constexpr void get_type_out(out_t & out, std::tuple<error_types_t...> const & in) const
    {
        constexpr auto count = meta::count<meta::list<error_types_t...>, error_type_t>::value;
        static_assert(count <= 1, "The same error type has been passed multiple times to max_error.");
        if constexpr (count > 0)
            std::get<error_type_t>(out) = std::get<error_type_t>(in);
    }

public:
    /*!\brief Adds to the configuration a max_error configuration element.
     * \relates seqan3::search_config_max_error
     * \param[in] cfg The configuration to be extended.
     * \param[in] ...error_types The maximum number of errors for each error type.
     * \returns A new configuration containing the max_error configuration element.
     */
    template <typename configuration_t, typename ... error_types_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    constexpr auto invoke(configuration_t && cfg, error_types_t && ...error_types) const
    {
        static_assert(is_valid_search_configuration_v<search_cfg::id::max_error, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(search_cfg::id::max_error));

        std::tuple in{error_types...};
        std::tuple<search_cfg::total<uint8_t>, search_cfg::substitution<uint8_t>,
                   search_cfg::insertion<uint8_t>, search_cfg::deletion<uint8_t>> out{0, 0, 0, 0};

        get_type_out<search_cfg::total<uint8_t>>(out, in);
        get_type_out<search_cfg::substitution<uint8_t>>(out, in);
        get_type_out<search_cfg::insertion<uint8_t>>(out, in);
        get_type_out<search_cfg::deletion<uint8_t>>(out, in);

        using error_types_list_t = meta::list<remove_cvref_t<error_types_t>...>;
        constexpr bool total_set = meta::in<error_types_list_t, search_cfg::total<uint8_t>>::value;
        constexpr bool other_error_types_set = meta::count<error_types_list_t, search_cfg::total<uint8_t>>::value !=
                                               meta::size<error_types_list_t>::value;

        // no specific error types specified
        if constexpr (!other_error_types_set)
        {
            // only total is set: set all to total
            if constexpr (total_set)
            {
                uint8_t const total_v{std::get<0>(out)};
                std::get<1>(out) = search_cfg::substitution<uint8_t>{total_v};
                std::get<2>(out) = search_cfg::insertion<uint8_t>{total_v};
                std::get<3>(out) = search_cfg::deletion<uint8_t>{total_v};
            }
        }
        // at least one specific error type specified
        else if constexpr (!total_set)
        {
            uint8_t const substitution_error = static_cast<uint8_t>(std::get<1>(out));
            uint8_t const insertion_error = static_cast<uint8_t>(std::get<2>(out));
            uint8_t const deletion_error = static_cast<uint8_t>(std::get<3>(out));
            // total not set. set it to sum of all error types
            std::get<0>(out) = search_cfg::total<uint8_t>{std::min<uint8_t>(255, substitution_error
                                                                               + insertion_error
                                                                               + deletion_error)};
        }

        search_config_max_error tmp{static_cast<std::tuple<uint8_t, uint8_t, uint8_t, uint8_t>>(out)};
        return std::forward<configuration_t>(cfg).push_front(std::move(tmp));
    }
};

//!\brief Helper template meta-function associated with detail::search_config_max_error.
//!\ingroup search_configuration
template <>
struct on_search_config<search_cfg::id::max_error>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename std::is_same<t, search_config_max_error>::type;
};

//!\brief Mapping from the detail::search_config_max_error type to it's corresponding seqan3::search_cfg::id.
//!\ingroup search_configuration
template <>
struct search_config_type_to_id<search_config_max_error>
{
    //!\brief The associated seqan3::search_cfg::id.
    static constexpr search_cfg::id value = search_cfg::id::max_error;
};
} // namespace seqan3::detail

namespace seqan3::search_cfg
{
/*!\brief A configuration element for the maximum number of errors across all error types (mismatches, insertions,
 *        deletions). This is an upper bound of errors independent from error numbers of specific error types.
 * \details An insertion corresponds to a base inserted into the query that does not occur in the text at the position,
 *          a deletion corresponds to a base deleted from the query sequence that does occur in the indexed text.
 *          Deletions at the beginning and at the end of the sequence are not considered during a search.
 * \ingroup search_configuration
 */
inline constexpr detail::search_config_max_error_adaptor max_error;

} // namespace seqan3::search_cfg
