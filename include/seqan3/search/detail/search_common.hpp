// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides data structures used by different search algorithms.
 */

#pragma once

#include <tuple>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

//!\brief Object grouping numbers of errors for different kind of error types.
struct search_param
{
    //!\brief Total number of errors (upper bound over all error types).
    uint8_t total{};
    //!\brief Total number of substitution errors.
    uint8_t substitution{};
    //!\brief Total number of insertion errors.
    uint8_t insertion{};
    //!\brief Total number of deletion errors.
    uint8_t deletion{};

    //!\brief Returns `true` if all member variables of `lhs` and `rhs` are equal, `false` otherwise.
    constexpr friend bool operator==(search_param const & lhs, search_param const & rhs) noexcept
    {
        return std::tie(lhs.total, lhs.substitution, lhs.insertion, lhs.deletion) ==
               std::tie(rhs.total, rhs.substitution, rhs.insertion, rhs.deletion);
    }

    //!\brief Returns `true` if any member variable of `lhs` and `rhs` are not equal, `false` otherwise.
    constexpr friend bool operator!=(search_param const & lhs, search_param const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
};

} // namespace seqan3::detail
