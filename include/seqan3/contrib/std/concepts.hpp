// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan-std/blob/main/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::stl concepts.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_CONCEPTS
#define SEQAN_STD_CONCEPTS

#include <ranges>

#if __cpp_lib_ranges >= 202110L
namespace seqan::stl::ranges
{

using std::ranges::viewable_range;

}
#else
#    include "detail/exposition_only.hpp"
namespace seqan::stl::ranges
{

template <class T>
concept viewable_range = std::ranges::range<T>
                      && ((std::ranges::view<std::remove_cvref_t<T>>
                           && std::constructible_from<std::remove_cvref_t<T>, T>)
                          || (!std::ranges::view<std::remove_cvref_t<T>>
                              && (std::is_lvalue_reference_v<T>
                                  || (std::movable<std::remove_reference_t<T>>
                                      && !seqan::stl::detail::is_initializer_list<std::remove_cvref_t<T>>))));

} // namespace seqan::stl::ranges

#endif

#endif // SEQAN_STD_CONCEPTS
