// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link core_configuration Core / Configuration submodule \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup core_configuration Configuration
 * \ingroup core
 * \see core
 * \brief Provides core functionality used to configure configurations.
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
 * \include test/snippet/core/configuration/configuration_setup.cpp
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
 * \include test/snippet/core/configuration/configuration_compatibility.cpp
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
 * \include test/snippet/core/configuration/configuration_combine.cpp

 * ### Access the data
 *
 * The configuration inherits from a std::tuple and exposes a tuple like interface using the standard
 * position-based and type-based `get` interfaces. The `get` interface was extended to also support
 * template template types as input template parameters to query the correct element:
 *
 * \include test/snippet/core/configuration/configuration_get.cpp
 *
 * The get interface returns a reference to the stored configuration element. In some cases, e.g. the implementor
 * of the actual algorithm, one wants to have an easy access to the actual value of the setting. Since, the
 * configuration must not contain all possible configuration elements the seqan3::configuration provides a
 * seqan3::configuration::get_or interface. Using this interface one can call get with a specific configuration element.
 * If this configuration element or a specialisation of it is already stored inside of the configuration, the respective
 * element is returned. Otherwise, the passed argument will be returned as the alternative.
 *
 * \include test/snippet/core/configuration/configuration_get_or.cpp
 */

#pragma once

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
