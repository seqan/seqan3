// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link utility utility module \endlink, which provides useful, stand-alone features
 *        that are not strongly coupled with SeqAn3.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup utility Utility
 * \brief Provides additional utility functionality used by multiple modules.
 *
 * \details
 * The utility module contains concepts, functions, traits and classes that
 * are independent of the remaining modules in SeqAn. These implementations
 * are considered external functionality, i.e. they could have been outsourced
 * into their own libraries.
 *
 * The utility module has no dependency to any other module except the \ref core
 * module.
 */

#pragma once

#include <seqan3/utility/char_operations/all.hpp>
#include <seqan3/utility/concept/all.hpp>
#include <seqan3/utility/container/all.hpp>
#include <seqan3/utility/math.hpp>
#include <seqan3/utility/parallel/all.hpp>
#include <seqan3/utility/range/all.hpp>
#include <seqan3/utility/simd/all.hpp>
#include <seqan3/utility/tuple/all.hpp>
#include <seqan3/utility/type_list/all.hpp>
#include <seqan3/utility/type_pack/all.hpp>
#include <seqan3/utility/type_traits/all.hpp>
#include <seqan3/utility/views/all.hpp>
