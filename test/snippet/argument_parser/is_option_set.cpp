// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"awesome-app", argc, argv}; // initialize

    int a{3};
    myparser.add_option(a, 'a', "awesome-parameter", "Please specify an integer.");

    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    if (myparser.is_option_set('a'))
        seqan3::debug_stream << "The user set option -a on the command line.\n";

    if (myparser.is_option_set("awesome-parameter"))
        seqan3::debug_stream << "The user set option --awesome-parameter on the command line.\n";

    // Asking for an option identifier that was not used before throws an error:
    // myparser.is_option_set("foo"); // throws seqan3::design_error

    return 0;
}
