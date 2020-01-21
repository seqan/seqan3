// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-Header for components of the algorithm submodule.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

/*!\defgroup algorithm Algorithm
 * \ingroup core
 * \brief Provides core functionality used to configure algorithms.
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
 * setting and must satisfy the seqan3::detail::config_element. The following snippet demonstrates a
 * basic setup for such configuration elements.
 *
 * \include test/snippet/core/algorithm/configuration_setup.cpp
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
 * \include test/snippet/core/algorithm/configuration_compatibility.cpp
 *
 * The type for the configuration element ids is used to overload the bool table. In the example above, both
 * elements can be combined with each other but not with themselves, to avoid inconsistent settings.
 *
 * ### Combining Configurations
 *
 * To enable easy creation of algorithm settings the seqan3::configuration supports a pipeable interface for the
 * different configuration elements which are added through the seqan3::pipeable_config_element base class.
 * The following snippet demonstrates how `bar` and `foo` can be combined.
 *
 * \include test/snippet/core/algorithm/configuration_combine.cpp

 * ### Access the data
 *
 * The configuration inherits from a std::tuple and exposes a tuple like interface using the standard
 * position-based and type-based `get` interfaces. The `get` interface was extended to also support
 * template template types as input template parameters to query the correct element:
 *
 * \include test/snippet/core/algorithm/configuration_get.cpp
 *
 * The get interface returns a reference to the stored configuration element. In some cases, e.g. the implementor
 * of the actual algorithm, one wants to have an easy access to the actual value of the setting. Since, the
 * configuration must not contain all possible configuration elements the seqan3::configuration provides a
 * seqan3::configuration::value_or interface, which provides direct access to the value of the respective
 * configuration element or uses a default value if the queried type is not contained in the configuration.
 *
 * \include test/snippet/core/algorithm/configuration_value_or.cpp
 */
 #include <seqan3/core/algorithm/bound.hpp>
#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/algorithm/configuration_utility.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
