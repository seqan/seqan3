// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

//!\author: Clemens Cords <clemens.cords@fu-berlin.de>

#include <benchmark/benchmark.h>
#include <cstring>
#include <fstream>

#include <seqan3/test/seqan2.hpp>
#if SEQAN3_HAS_SEQAN2
#include <seqan/seq_io.h>
#endif

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

constexpr unsigned default_seed = 1234u;

// ============================================================================
// generate fastq file
// ============================================================================
constexpr size_t default_sequence_length = 50; //length of nucleotide and quality sequence

std::string generate_fastq_file(size_t const entries_size)
{
    std::stringstream file{};
    std::string const id{"@name"};

    auto seed = default_seed;

    for (size_t i = 0; i < entries_size; ++i, ++seed)
    {
        auto seq = test::generate_sequence<dna5>(default_sequence_length, 0, seed);
        auto seq_as_chars = seq | views::to_char;

        auto qual = test::generate_sequence<phred42>(default_sequence_length, 0, seed);
        auto qual_as_chars = qual | views::to_char;

        file << id << '\n' << std::string(seq_as_chars.begin(), seq_as_chars.end()) << '\n' << '+' << '\n' << std::string(qual_as_chars.begin(), qual_as_chars.end()) << '\n';
    }

    return file.str();
}

// ============================================================================
// save file on disc temporarily
// ============================================================================
auto create_fastq_file(size_t const entries)
{
    // create temporary file, automatically removed on destruction
    test::tmp_filename fastq_file{"tmp.fastq"};
    auto fastq_file_path = fastq_file.get_path();

    // fill temporary file with a fastq file
    std::ofstream ostream{fastq_file_path};
    ostream << generate_fastq_file(entries);
    ostream.close();
    return fastq_file;
}

// ============================================================================
// write 3-line entry to stream as often as possible
// ============================================================================
void fastq_write_to_stream(benchmark::State & state)
{
    std::ostringstream ostream;
    sequence_file_output fout{ostream, format_fastq{}};

    std::string const id{"@name"};
    auto seq = test::generate_sequence<dna5>(default_sequence_length, 0, default_seed);
    auto qual = test::generate_sequence<phred42>(default_sequence_length, 0, default_seed);

    for (auto _ : state)
        fout.emplace_back(seq, id, qual);
}

// ============================================================================
// read dummy fastq file from a stream
// ============================================================================
void fastq_read_from_stream(benchmark::State & state)
{
    std::istringstream istream{generate_fastq_file(state.range(0))};
    sequence_file_input fin{istream, format_fastq{}};

    for (auto _ : state)
    {
        for (auto[seq, id, qual] : fin)
        {
            //silences [-Werror=unused-variable], no performance impact
            (void)seq, (void)id, (void)qual;
        }
    }


}

// ============================================================================
// seqan2 comparison
// ============================================================================
#if SEQAN3_HAS_SEQAN2
void fastq_read_from_stream_seqan2(benchmark::State & state)
{
    using namespace seqan;

    std::istringstream istream{};
    size_t entries_size = state.range(0);
    auto file = generate_fastq_file(entries_size);

    auto restart_iterator = [&istream, &file]()    // cf. format_fasta_benchmark
    {
        istream = std::istringstream{file};
        // same constant as seqan3 benchmark

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

        for (size_t i = 0; i < entries_size; ++i)
            readRecord(id, seq, qual, it, seqan::Fastq{});

        clear(id);
        clear(seq);
        clear(qual);
    }
}

#endif

// ============================================================================
// read dummy fastq file from temporary file on disk
// ============================================================================
void fastq_read_from_disk(benchmark::State & state)
{
    auto file = create_fastq_file(state.range(0));
    auto path = file.get_path();

    for (auto _ : state)
    {
        sequence_file_input fin{path};

        for (auto[seq, id, qual] : fin)
        {
            //silences [-Werror=unused-variable], no performance impact
            (void)seq, (void)id, (void)qual;
        }
    }
}

// ============================================================================
// seqan2 comparison
// ============================================================================
#if SEQAN3_HAS_SEQAN2

void fastq_read_from_disk_seqan2(benchmark::State & state)
{
    using namespace seqan;

    test::tmp_filename file_name{"tmp.fastq"};
    auto tmp_path = file_name.get_path();

    std::ofstream ostream{tmp_path};
    ostream << generate_fastq_file(state.range(0));
    ostream.close();

    StringSet<String<char>> ids{};
    StringSet<String<Dna5>> seqs{};
    StringSet<String<char>> quals{};

    for (auto _ : state)
    {
        SeqFileIn seqFileIn(tmp_path.c_str());
        readRecords(ids, seqs, quals, seqFileIn);

        clear(ids);
        clear(seqs);
        clear(quals);
    }
}

#endif

BENCHMARK(fastq_write_to_stream);

BENCHMARK(fastq_read_from_stream)->Arg(100);
BENCHMARK(fastq_read_from_disk)->Arg(100);

#if SEQAN3_HAS_SEQAN2
BENCHMARK(fastq_read_from_stream_seqan2)->Arg(100);
BENCHMARK(fastq_read_from_disk_seqan2)->Arg(100);
#endif

BENCHMARK(fastq_read_from_stream)->Arg(1000);
BENCHMARK(fastq_read_from_disk)->Arg(1000);

#if SEQAN3_HAS_SEQAN2
BENCHMARK(fastq_read_from_stream_seqan2)->Arg(1000);
BENCHMARK(fastq_read_from_disk_seqan2)->Arg(1000);
#endif

BENCHMARK(fastq_read_from_stream)->Arg(10000);
BENCHMARK(fastq_read_from_disk)->Arg(10000);

#if SEQAN3_HAS_SEQAN2
BENCHMARK(fastq_read_from_stream_seqan2)->Arg(100000);
BENCHMARK(fastq_read_from_disk_seqan2)->Arg(100000);
#endif

BENCHMARK_MAIN();
