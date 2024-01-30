// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, char const ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv}; // initialize

    //![validator_call]
    int myint;
    seqan3::value_list_validator my_validator{2, 4, 6, 8, 10};

    myparser.add_option(myint, 'i', "integer", "Give me a number.", seqan3::option_spec::standard, my_validator);
    //![validator_call]

    // an exception will be thrown if the user specifies an integer
    // that is not one of [2,4,6,8,10] (e.g. "./test_app -i 3")
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    seqan3::debug_stream << "integer given by user passed validation: " << myint << "\n";
    return 0;
}
