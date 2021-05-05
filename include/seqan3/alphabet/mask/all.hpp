// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the mask submodule; includes all headers from alphabet/mask/.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
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
 * Masks are useful in cases where an alphabet needs to be augmented with additional information.
 * A common use case is the introduction of don't-care positions (*masking*) in \link nucleotide nucleotide \endlink or
 * \link aminoacid aminoacid \endlink sequences without using ``'N'`` or ``'X'``, respectively.
 * Instead, the original alphabet is combined with seqan3::mask to create a seqan3::masked composite.
 * When printed via the seqan3::debug_stream, masked characters are displayed in lowercase.
 *
 * \include test/snippet/alphabet/mask/masked.cpp
 *
 * ### Types of masking
 *
 * There are two types of masking:
 *  * **hard-masking**, which converts masked characters to the ambiguous/unknown character (``'N'`` or ``'X'``)
 *  * **soft-masking**, which uses lowercase for masked characters
 *
 * ### Repeat masking
 *
 * The use of soft-masking was popularised by [RepeatMasker](https://www.repeatmasker.org/).
 * Interspersed repeats (transposons, retrotransposons and processed pseudogenes) and low complexity sequences are
 * denoted by lowercase characters. Note that larger repeats, such as segmental duplications, large tandem repeats and
 * whole gene duplications are generally not masked.
 */
