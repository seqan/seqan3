// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-include for the structure IO submodule.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

/*!\defgroup io_structure_file Structure File
 * \ingroup io
 * \brief Provides files and formats for handling structure data.
 *
 * \include{doc} doc/fragments/io_structure_input.md
 *
 * \include{doc} doc/fragments/io_structure_output.md
 */

#pragma once

#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/input_format_concept.hpp>
#include <seqan3/io/structure_file/input_options.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/io/structure_file/output_format_concept.hpp>
#include <seqan3/io/structure_file/output_options.hpp>
#include <seqan3/io/structure_file/record.hpp>
