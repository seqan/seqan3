// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \ingroup views
 * \brief [DEPRECATED] Auxiliary header for the \link views views submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \deprecated This header will be removed in 3.1.0; Please
 *             \#include seqan3/core/range/detail/adaptor_base.hpp, or
 *             \#include seqan3/core/range/detail/adaptor_from_functor.hpp, or
 *             \#include seqan3/core/range/detail/adaptor_for_view_without_args.hpp instead.
 */

#pragma once

#include <seqan3/core/range/detail/adaptor_base.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/range/detail/adaptor_for_view_without_args.hpp>

SEQAN3_DEPRECATED_HEADER(
    "This header is deprecated and will be removed in SeqAn-3.1.0; Please #include <seqan3/core/range/detail/{adaptor_base.hpp,adaptor_from_functor.hpp,adaptor_for_view_without_args.hpp}> instead.")
