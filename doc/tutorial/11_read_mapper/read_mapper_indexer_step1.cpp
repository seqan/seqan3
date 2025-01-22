// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

void run_program(std::filesystem::path const & reference_path, std::filesystem::path const & index_path)
{
    seqan3::debug_stream << "reference_file_path: " << reference_path << '\n';
    seqan3::debug_stream << "index_path           " << index_path << '\n';
}

struct cmd_arguments
{
    std::filesystem::path reference_path{};
    std::filesystem::path index_path{"out.index"};
};

void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "E. coli";
    parser.info.short_description = "Creates an index over a reference.";
    parser.info.version = "1.0.0";
    parser.add_option(args.reference_path,
                      'r',
                      "reference",
                      "The path to the reference.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta"}});
    parser.add_option(args.index_path,
                      'o',
                      "output",
                      "The output index file path.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"index"}});
}

//![main]
int main(int argc, char const ** argv)
{
    seqan3::argument_parser parser("Indexer", argc, argv);
    cmd_arguments args{};

    initialise_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }

    run_program(args.reference_path, args.index_path);

    return 0;
}
//![main]
