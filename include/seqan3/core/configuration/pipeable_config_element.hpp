// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::pipeable_config_element.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief Adds pipe interface to configuration elements.
 * \ingroup algorithm
 *
 * \details
 *
 * An empty base class which is used by the algorithm configurations to be identified as a configuration element and
 * to compose an algorithm configuration using the pipe-operator.
 */
struct pipeable_config_element
{};

} // namespace seqan3
