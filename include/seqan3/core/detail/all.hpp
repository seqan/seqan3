// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
