// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::structure_file_input_options.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \ingroup structure_file
 * \tparam seq_legal_alphabet The sequence legal alphabet exposed as type trait to the format.
 * \tparam structured_seq_combined Trait that exposes to the format whether seq and structure arguments are actually the
 * same/combined.
 */
template <typename seq_legal_alphabet, bool structured_seq_combined>
struct structure_file_input_options
{
    //!\brief Read the ID string only up until the first whitespace character.
    bool truncate_ids = false;
};

} // namespace seqan3
