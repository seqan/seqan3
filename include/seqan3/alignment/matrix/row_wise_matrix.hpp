// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::row_wise_matrix.
 */

#pragma once

#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>

namespace seqan3::detail
{

/*!\brief A matrix represented in a one-dimensional std::vector accessing the data in row-major-order.
 * \ingroup alignment_matrix
 * \implements seqan3::detail::matrix
 * \tparam value_t \copydoc seqan3::detail::two_dimensional_matrix::value_type
 */
template <typename value_t>
using row_wise_matrix = two_dimensional_matrix<value_t>;

} // namespace seqan3::detail
