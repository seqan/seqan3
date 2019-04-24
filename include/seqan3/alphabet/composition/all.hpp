// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Meta-header for the composition submodule; includes all headers from alphabet/composition/.
 */
#pragma once

#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/composition/union_composition.hpp>

/*!\defgroup composition Composition
 * \brief Provides data structures joining multiple alphabets into a single alphabet.
 * \ingroup alphabet
 *
 * \par Introduction
 *
 * Composition alphabets are special alphabets that allow you to combine existing alphabets into new ones. For example,
 * you can add new characters to existing alphabets by using seqan3::union_composition or combine alphabets with quality
 * information by using seqan3::cartesian_composition.
 *
 * We have currently two major composition alphabets:
 * * seqan3::cartesian_composition which roughly corresponds to the Cartesian product of the given types. It
 *   behaves similar to std::tuple, but it is specialised for alphabets.
 * * seqan3::union_composition which roughly corresponds to the Union of the given types. It behaves similar to
 *   std::variant, but it is specialised for alphabets.
 */
