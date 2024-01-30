// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/platform.hpp>
#if SEQAN3_WITH_CEREAL
//![complete]
#    include <fstream>

#    include <seqan3/argument_parser/all.hpp>
#    include <seqan3/core/debug_stream.hpp>
#    include <seqan3/io/sequence_file/input.hpp>
#    include <seqan3/search/all.hpp>
#    include <seqan3/search/fm_index/bi_fm_index.hpp>

#    include <cereal/archives/binary.hpp>

//! [solution]
struct reference_storage_t
{
    std::vector<std::string> ids;
    std::vector<std::vector<seqan3::dna5>> seqs;
};

void read_reference(std::filesystem::path const & reference_path, reference_storage_t & storage)
{
    seqan3::sequence_file_input reference_in{reference_path};
    for (auto && record : reference_in)
    {
        storage.ids.push_back(record.id());
        storage.seqs.push_back(record.sequence());
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
    // we need the alphabet and text layout before loading
    seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection> index;
    {
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }

    seqan3::sequence_file_input query_file_in{query_path};

    (void)sam_path; // Silence warning about unused parameter. We are not yet using `sam_path`.

    seqan3::configuration const search_config =
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{errors}}
        | seqan3::search_cfg::hit_all_best{};

    for (auto && record : query_file_in | std::views::take(20))
    {
        (void)storage; // Silence warning about unused parameter. We are not yet using `storage`.

        seqan3::debug_stream << "Hits:" << '\n';
        for (auto && result : search(record.sequence(), index, search_config))
            seqan3::debug_stream << result << '\n';
        seqan3::debug_stream << "======================" << '\n';
    }
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
//![complete]
#endif //SEQAN3_WITH_CEREAL
