// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Meta-header for the \link alignment_band band implementations \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \deprecated This header will be removed in 3.1.0; The contained functionality has been replaced by the
 *             seqan3::align_cfg::band_fixed_size configuration.
 */

#pragma once

/*!\defgroup alignment_band Band
 * \brief Data structures for computing banded sequence alignments.
 * \ingroup alignment
 * \see alignment
 *
 * \details
 *
 * SeqAn offers the computation of banded alignments to reduce the running time of the algorithm. This can be
 * helpful if the region in which the optimal alignment exists is known a priori. To specify the banded alignment
 * the developer can use the seqan3::align_cfg::band_fixed_size option.
 * This band configuration is initialised with a seqan3::align_cfg::lower_diagonal and a
 * seqan3::align_cfg::upper_diagonal. The upper diagonal must always be greater than or equal to the lower diagonal.
 * To choose the correct band parameters, imagine a matrix with the first sequence written on top and the second
 * sequence along the left vertical side. A negative value reflects a start of the diagonal within the vertical part,
 * while a positive value implies a start within the top part of this matrix at the respective position.
 *
 * \sa seqan3::align_cfg::band_fixed_size
 */

#include <seqan3/alignment/band/static_band.hpp>

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0. The contained functionality has been replaced by the seqan3::align_cfg::band_fixed_size configuration.")
