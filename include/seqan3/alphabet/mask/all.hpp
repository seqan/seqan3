// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Meta-header for the mask submodule; includes all headers from alphabet/mask/.
 */

#pragma once

#include <seqan3/alphabet/mask/mask.hpp>
#include <seqan3/alphabet/mask/masked.hpp>

/*!\defgroup mask Mask
 * \brief Provides the mask alphabet and functionality for creating masked composites.
 * \ingroup alphabet
 *
 * \details
 *
 * ### Introduction
 *
 * Masks are useful as tuple composites when one wants to create a masked alphabet with
 * don't care positions, but does not want to use the seqan3::dna15 **N** or
 * seqan3::aa27 **X** because of loss of information. It will instead mark the specified characters as masked,
 * and display them as lowercase representations when printed.\n
 * There are two types of masking: "hard-masking" which converts to the UNKNOWN character and
 * "soft-masking", which is visualised by using lower-case instead of upper-case.
 * However because regular nucleotide and aminoacid alphabets discard case on assignment,
 * one needs to create additional alphabets to preserve this information (if desired).\n
 * This alphabet in itself is not useful to users directly, but instead the composite seqan3::masked may be used to
 * transform another alphabet into a new alphabet that can represent the original alphabet plus masking information.
 */
