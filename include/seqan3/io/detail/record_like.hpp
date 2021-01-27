// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::record_like.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/record.hpp>
#include <seqan3/utility/tuple/concept.hpp>

namespace seqan3::detail
{

//!\brief Helper struct to implement seqan3::detail::record_like
template <typename record_t>
struct is_derived_from_record
{
private:
    //!\overload
    static std::false_type derived_from(...);

    //!\brief Helper function to determine whether the given record_t is derived_from seqan3::record.
    template <typename ...args_t>
    static std::true_type derived_from(seqan3::record<args_t...> &);
public:

    //!\brief Whether the given record_t is derived_from seqan3::record.
    static constexpr bool value = decltype(derived_from(std::declval<record_t &>())){};
};

/*!\interface seqan3::detail::record_like <>
 * \brief The concept for a type that models a record.
 * \ingroup io
 */
//!\cond
template <typename record_t>
SEQAN3_CONCEPT record_like = tuple_like<record_t> &&
                             is_derived_from_record<std::remove_cvref_t<record_t>>::value;
//!\endcond
} // namespace seqan3::detail
