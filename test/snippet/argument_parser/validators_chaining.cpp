// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>

#include <seqan3/argument_parser/all.hpp>

int main(int argc, char const ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv}; // initialize

    //![validator_call]
    std::string file_name;

    seqan3::regex_validator absolute_path_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"};
    seqan3::input_file_validator my_file_ext_validator{{"sa", "so"}};

    myparser.add_option(file_name,
                        'f',
                        "file",
                        "Give me a file name with an absolute path.",
                        seqan3::option_spec::standard,
                        absolute_path_validator | my_file_ext_validator);
    //![validator_call]

    // an exception will be thrown if the user specifies a file name
    // that is not an absolute path or does not have one of the file extension [sa,so]
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    std::cout << "filename given by user passed validation: " << file_name << "\n";
    return 0;
}
