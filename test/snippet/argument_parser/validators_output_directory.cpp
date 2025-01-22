// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, char const ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv}; // initialize

    //! [validator_call]
    std::filesystem::path mydir{};

    myparser.add_option(mydir,
                        'd',
                        "dir",
                        "The output directory for storing the files.",
                        seqan3::option_spec::standard,
                        seqan3::output_directory_validator{});
    //! [validator_call]

    // an exception will be thrown if the user specifies a directory that cannot be created by the filesystem either
    // because the parent path does not exists or the path has insufficient write permissions.
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    seqan3::debug_stream << "directory given by user passed validation: " << mydir << "\n";
    return 0;
}
