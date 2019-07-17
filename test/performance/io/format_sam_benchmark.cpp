// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

// !\author Clemens Cords (clemens.cords@fu-berlin.de)

#include <cctype>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <streambuf>

#include <benchmark/benchmark.h>

#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/test/tmp_filename.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/bam_io.h>
#endif

using namespace seqan3;

constexpr auto default_seed = 1234u;
size_t const reference_length = 50;

// ============================================================================
// generate mutated sequence from reference
// ============================================================================
enum class Mutation_type : uint8_t
{
    DELETION,
    INSERTION,
    SUBSTITUTION
};

float const mutation_probability = 0.1;

template<Alphabet alphabet_t>
std::vector<alphabet_t> mutate_sequence(std::vector<alphabet_t> input_sequence,
                                        size_t seed = default_seed,
                                        float mutation_probability = mutation_probability)
{
    std::mt19937 rng_engine(seed);

    // c.f. for loop below
    std::uniform_real_distribution<float> cointoss_distribution(0, 1);
    std::uniform_int_distribution<uint8_t> mutation_distribution(0, 3 -1); // 3 = n values in Mutation_type;
    std::uniform_int_distribution<uint8_t> base_distribution(0, alphabet_size<alphabet_t> -1);

    std::vector<alphabet_t> mutated_sequence{};

    for (std::size_t i = 0; i < input_sequence.size(); ++i)
    {
        // cointoss if letter should be mutated at all
        if (cointoss_distribution(rng_engine) > mutation_probability)
        {
            mutated_sequence.push_back(input_sequence[i]);
        }
        else
        {
            // determine type of mutation and mutate
            Mutation_type mutation_type = static_cast<Mutation_type>(mutation_distribution(rng_engine));

            switch (mutation_type)
            {
                case Mutation_type::DELETION:
                    break; // skip

                case Mutation_type::INSERTION :
                    mutated_sequence.push_back(alphabet_t{}.assign_rank(mutation_distribution(rng_engine)));
                    mutated_sequence.push_back(input_sequence[i]);
                    break;

                case Mutation_type::SUBSTITUTION:
                    mutated_sequence.push_back(alphabet_t{}.assign_rank(mutation_distribution(rng_engine)));
                    break;
            }
        }
    }

    return mutated_sequence;
}

// ============================================================================
// generate sam file from randomly generated sequence and store in file_dict
// ============================================================================
std::map<std::size_t, std::string> file_dict;
std::string const default_query_id = "query_";
std::string reference_id = "reference_id";

static std::string get_file (std::size_t n_queries)
{
    if (file_dict.find(n_queries) == file_dict.end())
    {
        std::string file{};

        // generate sequences
        auto reference = test::generate_sequence<dna4>(reference_length, 0, default_seed);
        std::vector<std::vector<dna4>> queries;

        size_t seed = default_seed;
        for (size_t i = 0; i < n_queries; ++i, ++seed)
        {
            auto query = mutate_sequence<dna4>(reference, seed);
            queries.push_back(query);
        }

        // align
        auto config = align_cfg::mode{global_alignment} |
                      align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-2}}} |
                      align_cfg::gap{gap_scheme{gap_score{-4}}} |
                      align_cfg::aligned_ends{free_ends_all} |
                      align_cfg::result{with_alignment};

        std::stringstream sam_stream;
        format_sam format{};

        alignment_file_output sam_out{sam_stream,
                                      format,
                                      fields< field::SEQ,
                                              field::ID,
                                              field::REF_ID,
                                              field::REF_SEQ,
                                              field::ALIGNMENT,
                                              field::MAPQ>{}};

        size_t i = 0;
        for (auto query : queries)
        {
            auto & result = *(align_pairwise(std::tie(reference, query), config).begin());
            auto aligned_sequence = result.alignment();
            size_t map_qual = 60u + result.score();
            std::string current_query_id = default_query_id + std::to_string(i);
            ++i;
            sam_out.emplace_back(query, current_query_id, reference_id, reference, aligned_sequence, map_qual);

        }

        // convert to string
        file = sam_stream.str();

        file_dict[n_queries] = file;
    }

    return file_dict[n_queries];
}

// ============================================================================
// read sam files directly from filestream
// ============================================================================
void sam_file_read_from_stream(benchmark::State &state)
{
    size_t n_queries = state.range(0);
    std::string sam_file = get_file(n_queries);

    std::istringstream istream{sam_file};

    for (auto _ : state)
    {
        // refill stream (same constant for seqan2 benchmark)
        istream.clear();
        istream.seekg(0, std::ios::beg);

        alignment_file_input fin{istream, format_sam{}};

        // read all records and store in internal buffer
        auto it = fin.begin();
        while (it != fin.end())
            ++it;
    }
}

BENCHMARK(sam_file_read_from_stream)->Arg(50);
BENCHMARK(sam_file_read_from_stream)->Arg(500);
BENCHMARK(sam_file_read_from_stream)->Arg(5000);

#if SEQAN3_HAS_SEQAN2

// ============================================================================
// seqan2 read from stream
// ============================================================================
void seqan2_sam_file_read_from_stream(benchmark::State &state)
{
    using namespace seqan;

    size_t n_queries = state.range(0);
    std::string sam_file = get_file(n_queries);

    // create temporary BamFileIn and read from disk
    // to get the stupid context...
    test::tmp_filename file_name{"tmp.sam"};
    {
        std::ofstream ostream{file_name.get_path()};
        ostream << get_file(n_queries);
        ostream.close();
    }

    BamHeader tmp_header;
    BamFileIn tmp_bam_file_in(file_name.get_path().c_str());
    readHeader(tmp_header, tmp_bam_file_in);
    auto cxt = context(tmp_bam_file_in);

    BamAlignmentRecord record;
    BamHeader header;

    std::istringstream istream{sam_file};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);

        auto it = seqan::Iter<std::istringstream, seqan::StreamIterator<seqan::Input> >(istream);

        readHeader(header, cxt, it, Sam());

        for (size_t i = 0; i < n_queries; ++i)
        {
            readRecord(record, cxt, it, Sam());
            clear(record);
        }

        clear(header);
    }
}

BENCHMARK(seqan2_sam_file_read_from_stream)->Arg(50);
BENCHMARK(seqan2_sam_file_read_from_stream)->Arg(500);
BENCHMARK(seqan2_sam_file_read_from_stream)->Arg(5000);

#endif // has seqan2

// ============================================================================
// seqan3 read from disk
// ============================================================================
void sam_file_read_from_disk(benchmark::State &state)
{
    // create temporary file
    test::tmp_filename file_name{"tmp.sam"};
    auto tmp_path = file_name.get_path();

    size_t n_queries = state.range(0);

    std::ofstream ostream{tmp_path};
    ostream << get_file(n_queries);
    ostream.close();

    for (auto _ : state)
    {
        alignment_file_input fin{tmp_path};

        // read all records and store in internal buffer
        auto it = fin.begin();
        while (it != fin.end())
            ++it;
    }
}

BENCHMARK(sam_file_read_from_disk)->Arg(50);
BENCHMARK(sam_file_read_from_disk)->Arg(500);
BENCHMARK(sam_file_read_from_disk)->Arg(5000);

#if SEQAN3_HAS_SEQAN2

// ============================================================================
// seqan2 read from disk
// ============================================================================
void seqan2_sam_file_read_from_disk(benchmark::State &state)
{
    using namespace seqan;

    // tmp file
    test::tmp_filename file_name{"tmp.sam"};
    auto tmp_path = file_name.get_path();

    size_t n_queries = state.range(0);

    std::ofstream ostream{tmp_path};
    ostream << get_file(n_queries);
    ostream.close();

    BamHeader header;
    BamAlignmentRecord record;

    for (auto _ : state)
    {
        BamFileIn bamFileIn(tmp_path.c_str());

        readHeader(header, bamFileIn);

        while (!atEnd(bamFileIn))
            readRecord(record, bamFileIn);

        clear(header);
        clear(record);
    }
}

BENCHMARK(seqan2_sam_file_read_from_disk)->Arg(50);
BENCHMARK(seqan2_sam_file_read_from_disk)->Arg(500);
BENCHMARK(seqan2_sam_file_read_from_disk)->Arg(5000);

#endif // has seqan2

BENCHMARK_MAIN();
