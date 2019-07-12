#include <seqan3/core/platform.hpp>
#if SEQAN3_WITH_CEREAL
//![complete]
#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/algorithm/all.hpp>

using namespace seqan3;

//! [solution]
struct reference_storage_t
{
    std::vector<std::string> ids;
    std::vector<std::vector<dna5>> seqs;
};

void read_reference(std::filesystem::path const & reference_path,
                    reference_storage_t & storage)
{
    sequence_file_input reference_in{reference_path};
    for (auto & [seq, id, qual] : reference_in)
    {
        storage.ids.push_back(std::move(id));
        storage.seqs.push_back(std::move(seq));
    }
}

//! [map_reads]
void map_reads(std::filesystem::path const & query_path,
               std::filesystem::path const & index_path,
               std::filesystem::path const & sam_path,
               reference_storage_t & storage,
               uint8_t const errors)
//! [map_reads]
{
    bi_fm_index<text_layout::collection> index; // we need to know if we work on a text collection before loading
    {
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }

    sequence_file_input query_in{query_path};

    configuration const search_config = search_cfg::max_error{search_cfg::total{errors}} |
                                        search_cfg::mode{search_cfg::all_best};

    for (auto & [query, id, qual] : query_in | view::take(20))
    {
        auto positions = search(query, index, search_config);
        debug_stream << "id:           " << id << '\n';
        debug_stream << "positions:    " << positions << '\n';
        debug_stream << "======================" << '\n';
    }
    (void) sam_path; // prevent unused parameter warning
    (void) storage; // prevent unused parameter warning
}

void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & index_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
    reference_storage_t storage{};
    read_reference(reference_path, storage);
    map_reads(query_path, index_path, sam_path, storage, errors);
}
//! [solution]

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
//![complete]
#endif //SEQAN3_WITH_CEREAL
