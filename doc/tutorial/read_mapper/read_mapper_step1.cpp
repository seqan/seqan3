#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & index_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
    debug_stream << "reference_path: " << reference_path << '\n';
    debug_stream << "query_path:     " << query_path << '\n';
    debug_stream << "index_path      " << index_path << '\n';
    debug_stream << "sam_path        " << sam_path << '\n';
    debug_stream << "errors:         " << errors << '\n';
}

struct cmd_arguments
{
    std::filesystem::path reference_path{};
    std::filesystem::path query_path{};
    std::filesystem::path index_path{};
    std::filesystem::path sam_path{"out.sam"};
    uint8_t errors{0};
};

void initialise_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "E. coli";
    parser.info.short_description = "Map reads against a reference.";
    parser.info.version = "1.0.0";
    parser.add_option(args.reference_path, 'r', "reference", "The path to the reference.", option_spec::REQUIRED,
                      input_file_validator{{"fa","fasta"}});
    parser.add_option(args.query_path, 'q', "query", "The path to the query.", option_spec::REQUIRED,
                      input_file_validator{{"fq","fastq"}});
    parser.add_option(args.index_path, 'i', "index", "The path to the index.", option_spec::REQUIRED,
                      input_file_validator{{"index"}});
    parser.add_option(args.sam_path, 'o', "output", "The output SAM file path.", option_spec::DEFAULT,
                      output_file_validator{{"sam"}});
    parser.add_option(args.errors, 'e', "error", "Maximum allowed errors.", option_spec::DEFAULT,
                      arithmetic_range_validator{0, 4});
}

//![main]
int main(int argc, char const ** argv)
{
    argument_parser parser("Mapper", argc, argv);
    cmd_arguments args{};

    initialise_argument_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }

    run_program(args.reference_path, args.query_path, args.index_path, args.sam_path, args.errors);

    return 0;
}
//![main]
