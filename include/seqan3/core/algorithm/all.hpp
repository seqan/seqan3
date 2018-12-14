// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Meta-Header for components of the algorithm submodule.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

/*!\defgroup algorithm Algorithm
 * \ingroup core
 * \brief Contains core functionality used to configure algorithms.
 *
 * \details
 *
 * In SeqAn there are many algorithms, e.g. alignment or search algorithms, that can be configured through
 * many different settings and policies that alter the execution of the respective algorithm.
 * These configurations can be orthogonal or might be mutually exclusive and can make interfaces very difficult to use.
 * This module provides a basic system to manage the configurations of algorithms using a unified interface.
 *
 * ### Usage
 *
 * The basis of any algorithm configuration are configuration elements. These are objects that handle a specific
 * setting and must satisfy the seqan3::detail::config_element_concept. The following snippet demonstrates a
 * basic setup for such configuration elements.
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp configuration_setup
 *
 * Here, two classes with the name `bar` and `foo` are created. They can be normal or template classes and must
 * inherit from seqan3::pipeable_config_element and contain a member with the name `value` to satisfy the respective
 * concept. The separate `value` member is used for a proper encapsulation from the actual setting parameter.
 * For example the alignment algorithms require a scoring scheme, but the scoring scheme itself should not be pipeable
 * with other settings.
 *
 * In addition an enum type was defined that will be used later to allow for compatibility checks when combining
 * settings in a configuration object. This enum assigns every configuration element a unique id.
 * To provide compatibility checks the `seqan3::detail::compatibility_table` must be overloaded for the specific
 * algorithm configuration.
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp compatibility
 *
 * The type for the configuration element ids is used to overload the bool table. In the example above, both
 * elements can be combined with each other but not with themselves, to avoid inconsistent settings.
 *
 * ### Combining Configurations
 *
 * To enable easy creation of algorithm settings the seqan3::configuration supports a pipeable interface for the
 * different configuration elements which is added through the seqan3::pipeable_config_element base class.
 * The following snippet demonstrates how `bar` and `foo` can be combined.
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp combine

 * ### Access the data
 *
 * The configuration inherits from a std::tuple and exposes a tuple like interface using the standard
 * position-based and type-based `get` interfaces. The `get` interface was extended to also support
 * template template types as input template parameters to query the correct element:
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp get
 *
 * The get interface returns a reference to the stored configuration element. In some cases, e.g. the implementor
 * of the actual algorithm, one wants to have an easy access to the actual value of the setting. Since, the
 * configuration must not contain all possible configuration elements the seqan3::configuration provides a
 * seqan3::configuration::value_or interface, which provides direct access to the value of the respective
 * configuration element or uses a default value if the queried type is not contained in the configuration.
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp value_or
 */
#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/algorithm/configuration_utility.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
