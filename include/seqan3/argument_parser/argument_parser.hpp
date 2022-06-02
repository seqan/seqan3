// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::argument_parser class.
 */

#pragma once

#include <seqan3/argument_parser/auxiliary.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/core/detail/test_accessor.hpp>

#include <sharg/argument_parser.hpp>

namespace seqan3
{

using namespace sharg;

/*!\brief The SeqAn command line parser.
 * \ingroup argument_parser
 */
class argument_parser /* \cond */ : public sharg::argument_parser /* \endcond */
{
    //!\cond
    // the documentation is included by the sharg parser.
    // Only API that is not present in sharg needs to be documented.
    using base_t = sharg::argument_parser;

    friend struct ::seqan3::detail::test_accessor;

public:
    argument_parser() = delete;
    argument_parser(argument_parser const &) = default;
    argument_parser & operator=(argument_parser const &) = default;
    argument_parser(argument_parser &&) = default;
    argument_parser & operator=(argument_parser &&) = default;

    argument_parser(std::string const app_name,
                    int const argc,
                    char const * const * const argv,
                    update_notifications version_updates = update_notifications::on,
                    std::vector<std::string> subcommands = {}) :
        base_t{app_name, argc, argv, version_updates, subcommands}
    {}

    // needs to be overloaded s.t. the correct type is returned (not the sharg::argument_parser)
    argument_parser & get_sub_parser()
    {
        return static_cast<argument_parser &>(base_t::get_sub_parser());
    }
    //!\endcond
};

} // namespace seqan3
