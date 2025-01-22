// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & index_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
    seqan3::debug_stream << "reference_path: " << reference_path << '\n';
    seqan3::debug_stream << "query_path:     " << query_path << '\n';
    seqan3::debug_stream << "index_path      " << index_path << '\n';
    seqan3::debug_stream << "sam_path        " << sam_path << '\n';
    seqan3::debug_stream << "errors:         " << errors << '\n';
}

struct cmd_arguments
{
    std::filesystem::path reference_path{};
    std::filesystem::path query_path{};
    std::filesystem::path index_path{};
    std::filesystem::path sam_path{"out.sam"};
    uint8_t errors{0};
};

void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "E. coli";
    parser.info.short_description = "Map reads against a reference.";
    parser.info.version = "1.0.0";
    parser.add_option(args.reference_path,
                      'r',
                      "reference",
                      "The path to the reference.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta"}});
    parser.add_option(args.query_path,
                      'q',
                      "query",
                      "The path to the query.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fq", "fastq"}});
    parser.add_option(args.index_path,
                      'i',
                      "index",
                      "The path to the index.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"index"}});
    parser.add_option(args.sam_path,
                      'o',
                      "output",
                      "The output SAM file path.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"sam"}});
    parser.add_option(args.errors,
                      'e',
                      "error",
                      "Maximum allowed errors.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0, 4});
}

//![main]
int main(int argc, char const ** argv)
{
    seqan3::argument_parser parser("Mapper", argc, argv);
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

    run_program(args.reference_path, args.query_path, args.index_path, args.sam_path, args.errors);

    return 0;
}
//![main]
