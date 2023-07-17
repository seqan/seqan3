// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for implementation details in the core module.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <bit>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/core/detail/deferred_crtp_base.hpp>
#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
