#include <seqan3/core/platform.hpp>
#if SEQAN3_WITH_CEREAL
//![complete]
#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

struct reference_storage_t
{
    std::vector<std::string> ids;
    std::vector<std::vector<seqan3::dna5>> seqs;
};

//! [solution]
void read_reference(std::filesystem::path const & reference_path,
                    reference_storage_t & storage)
{
    seqan3::sequence_file_input reference_in{reference_path};
    for (auto & [seq, id, qual] : reference_in)
    {
        storage.ids.push_back(std::move(id));
        storage.seqs.push_back(std::move(seq));
    }
}


//! [create_index]
void create_index(std::filesystem::path const & index_path,
                  reference_storage_t & storage)
//! [create_index]
{
    seqan3::bi_fm_index index{storage.seqs};
    {
        std::ofstream os{index_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }
}

void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & index_path)
{
    reference_storage_t storage{};
    read_reference(reference_path, storage);
    create_index(index_path, storage);
}
//! [solution]

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
    parser.add_option(args.reference_path, 'r', "reference", "The path to the reference.",
                      seqan3::option_spec::REQUIRED,
                      seqan3::input_file_validator{{"fa","fasta"}});
    parser.add_option(args.index_path, 'o', "output", "The output index file path.",
                      seqan3::option_spec::DEFAULT,
                      seqan3::output_file_validator{{"index"}});
}

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
//![complete]
#endif //SEQAN3_WITH_CEREAL
