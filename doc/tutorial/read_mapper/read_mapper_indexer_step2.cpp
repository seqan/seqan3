#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

//! [solution]
//! [reference_storage_t]
struct reference_storage_t
{
    std::vector<std::string> ids;
    std::vector<std::vector<dna5>> seqs;
};
//! [reference_storage_t]

//! [read_reference]
void read_reference(std::filesystem::path const & reference_path,
                    reference_storage_t & storage)
//! [read_reference]
{
    sequence_file_input reference_in{reference_path};
    for (auto & [seq, id, qual] : reference_in)
    {
        storage.ids.push_back(std::move(id));
        storage.seqs.push_back(std::move(seq));
    }
    debug_stream << "Reference IDs: " << storage.ids << '\n';
}

void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & index_path)
{
    debug_stream << "reference_file_path: " << reference_path << '\n';
    debug_stream << "index_path           " << index_path << '\n';
    reference_storage_t storage{};
    read_reference(reference_path, storage);
}
//! [solution]

struct cmd_arguments
{
    std::filesystem::path reference_path{};
    std::filesystem::path index_path{"out.index"};
};

void initialise_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "E. coli";
    parser.info.short_description = "Creates an index over a reference.";
    parser.info.version = "1.0.0";
    parser.add_option(args.reference_path, 'r', "reference", "The path to the reference.", option_spec::REQUIRED,
                      input_file_validator{{"fa","fasta"}});
    parser.add_option(args.index_path, 'o', "output", "The output index file path.", option_spec::DEFAULT,
                      output_file_validator{{"index"}});
}

int main(int argc, char const ** argv)
{
    argument_parser parser("Indexer", argc, argv);
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

    run_program(args.reference_path, args.index_path);

    return 0;
}
