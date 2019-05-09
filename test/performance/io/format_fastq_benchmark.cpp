// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

// written for the purpose of the Bch Softwarepraktikum, early 2019 FU Berlin
// @author: Clemens Cords <clemenscords@fu-berlin.de>

#include <benchmark/benchmark.h>
#include <cstring>
#include <cmath>
#include <fstream>

#if __has_include(<seqan/seq_io.h>)
    #include <seqan/seq_io.h>
#endif

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/std/filesystem>

using namespace seqan3;

unsigned int const SEED = 1234;

// ============================================================================
// generate fastq file
// ============================================================================
static constexpr std::size_t SEQUENCE_LENGTH = 50; // length of nucleotide and quality sequence
static constexpr std::size_t N_ENTRIES_IN_FILE = 4069; // number of 3-line entries

// run-time execution workaround because test::generate_sequence does not work at compile time
std::string fastq_file;
bool has_been_init = false;

std::string get_file()
{
    if (!has_been_init) {
        std::string file{};
        std::string const id{"@name"};

        for (size_t i = 0; i < N_ENTRIES_IN_FILE; ++i)
        {
            auto seq = test::generate_sequence<dna5>(SEQUENCE_LENGTH, 0, SEED);
            std::string seq_string = seq | view::to_char;

            auto quality = test::generate_sequence<phred42>(SEQUENCE_LENGTH, 0, SEED);
            std::string quality_string = quality | view::to_char;
            file += id + '\n' + seq_string + '\n' + '+' + '\n' + quality_string + '\n';
        }

        has_been_init = true;
        fastq_file = file;
    }

    return fastq_file;
}

// ============================================================================
// try to write 3-line entry to stream as often as possible
// ============================================================================
void fastq_write(benchmark::State & state)
{
    std::ostringstream ostream;

    sequence_file_format_fastq format;
    sequence_file_output_options options{};

    std::string const id{"@name"};
    auto seq = test::generate_sequence<dna5>(SEQUENCE_LENGTH, 0, SEED);
    auto quality = test::generate_sequence<phred42>(SEQUENCE_LENGTH, 0, SEED);

    for (auto _ : state)
    {
        format.write(ostream, options, seq, id, quality);
    }
}

BENCHMARK(fastq_write);


// ============================================================================
// read dummy fastq file from a stream
// ============================================================================
void fastq_read_from_stream(benchmark::State & state)
{
    sequence_file_format_fastq format;
    sequence_file_input_options<dna5, false> const options{};

    std::string id{};
    std::vector<dna5> seq{};
    std::vector<phred42> qual{};

    for (auto _ : state)
    {
        std::istringstream istream{get_file()};
        // refilling stream skews benchmark but unavoidable

        format.read(istream, options, seq, id, qual);

        id.clear();
        seq.clear();
        qual.clear();
    }
}

BENCHMARK(fastq_read_from_stream);

// ============================================================================
// seqan2 comparison
// ============================================================================
#if __has_include(<seqan/seq_io.h>)

void fastq_read_from_stream_seqan2(benchmark::State & state)
{
    using namespace seqan;

    std::istringstream istream{};

    auto restart_iterator = [&istream]()    // c.f. format_fasta_benchmark
    {
        istream = std::istringstream{get_file()};
        seqan::VirtualStream<char, Input> comp{};
        open(comp, istream);
        return directionIterator(comp, Input());
    };

    String<char> id{};
    String<Dna5> seq{};
    String<char> qual{};

    for (auto _ : state)
    {
        auto it = restart_iterator();
        // should be same overhead as istream(file) in equivalent seqan3 benchmark
        // thus skews benchmark but still makes for a valid comparison between seqan2 and 3

        readRecord(id, seq, qual, it, seqan::Fastq{});
        clear(id);
        clear(seq);
        clear(qual);
    }
}

BENCHMARK(fastq_read_from_stream_seqan2);
#endif

// ============================================================================
// read dummy fastq file from temporary file on disc
// ============================================================================
void fastq_read_from_disc(benchmark::State & state)
{
    sequence_file_format_fastq format;

    // create temporary file on disc
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();
    std::string const file_name = tmp_dir/"tmp.fastq";

    std::ofstream ostream{file_name};
    bool open_success = ostream.is_open();

    if (!open_success)
    {
        std::perror("Error creating temporary file \"tmp.fastq\" for fastq_read_from_disc benchmark.");
        std::cout << "aborting..." << std::endl;
        exit(1);
    }

    ostream << get_file();
    ostream.close();

    // benchmark
    std::string id{};
    std::vector<dna5> seq{};
    std::vector<phred42> qual{};

    for (auto _ : state)
    {
        sequence_file_input fin{file_name};

        // read all records and store in internal buffer
        auto it = fin.begin();
        while (it != fin.end())
            ++it;
    }

    // remove temporary file
    const int remove_success = std::remove(file_name.c_str());
    if (remove_success != 0)    //sic
    {
        std::perror("Error removing temporary file \"tmp.fastq\" for fastq_read_from_disc benchmark");
        std::cout << "aborting..." << std::endl;
        exit(1);
    }
}

BENCHMARK(fastq_read_from_disc);

#if __has_include(<seqan/seq_io.h>)
void fastq_read_from_disc_seqan2(benchmark::State & state)
{
    using namespace seqan;

    // create
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();
    std::string const file_name = tmp_dir/"tmp.fastq";

    std::ofstream ostream{file_name};
    bool open_success = ostream.is_open();

    if (!open_success)
    {
        std::perror("Error creating temporary file \"tmp.fastq\" for fastq_read_from_disc_seqan2.");
        std::cout << "aborting..." << std::endl;
        exit(1);
    }

    ostream << get_file();
    ostream.close();

    // benchmark
    StringSet<String<char>> ids{};
    StringSet<String<Dna5>> seqs{};
    StringSet<String<char>> quals{};

    for (auto _ : state)
    {
        SeqFileIn seqFileIn(file_name.c_str());
        readRecords(ids, seqs, quals, seqFileIn);

        clear(ids);
        clear(seqs);
        clear(quals);
    }

    // delete
    const int remove_success = std::remove(file_name.c_str());
    if (remove_success != 0)
    {
        std::perror("Error removing temporary file \"tmp.fastq\" for fastq_read_from_disc_seqan2.");
        std::cout << "aborting..." << std::endl;
        exit(1);
    }
}

BENCHMARK(fastq_read_from_disc_seqan2);

#endif

/*
// ============================================================================
// read dummy fastq file ignoring all parameters
// ============================================================================
void fastq_read_ignore_everything(benchmark::State & state)
{
    auto const fastq_file = generate_dummy_fastq_file(POWER);

    sequence_file_format_fastq format;
    sequence_file_input_options<dna5, false> const options{};

    for (auto _ : state)
    {
        state.PauseTiming();
        std::istringstream istream{fastq_file};
        state.ResumeTiming();

        for (size_t i = 0; i < iterations_per_run; ++i)
            format.read(istream, options, std::ignore, std::ignore, std::ignore);
    }
}

BENCHMARK(fastq_read_ignore_everything);
*/

BENCHMARK_MAIN();
