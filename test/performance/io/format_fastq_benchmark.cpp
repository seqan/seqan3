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
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/view/to_char.hpp>
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
    std::string file{};
    std::string const id{"@name"};

    auto seed = default_seed;

    for (size_t i = 0; i < entries_size; ++i, ++seed)
    {
        auto seq = test::generate_sequence<dna5>(default_sequence_length, 0, seed);
        std::string seq_string = seq | view::to_char;

        auto qual = test::generate_sequence<phred42>(default_sequence_length, 0, seed);
        std::string qual_string = qual | view::to_char;

        file += id + '\n' + seq_string + '\n' + '+' + '\n' + qual_string + '\n';
    }

    return file;
}

// ============================================================================
// write 3-line entry to stream as often as possible
// ============================================================================
void fastq_write_to_stream(benchmark::State & state)
{
    std::ostringstream ostream;

    sequence_file_format_fastq format;
    sequence_file_output_options options{};

    std::string const id{"@name"};
    auto seq = test::generate_sequence<dna5>(default_sequence_length, 0, default_seed);
    auto qual = test::generate_sequence<phred42>(default_sequence_length, 0, default_seed);

    for (auto _ : state)
    {
        format.write(ostream, options, seq, id, qual);
    }
}

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

    size_t n_entries_in_file = state.range(0);
    auto file = generate_fastq_file(n_entries_in_file);

    for (auto _ : state)
    {
        std::istringstream istream{file};
        // refilling streams introduces the same constant to
        // both seqan3 and 2 benchmark runtime so they are still comparable

        for (size_t i = 0; i < n_entries_in_file; ++i)
            format.read(istream, options, seq, id, qual);

        id.clear();
        seq.clear();
        qual.clear();
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
    size_t n_entries_in_file = state.range(0);
    auto file = generate_fastq_file(n_entries_in_file);

    auto restart_iterator = [&istream, &file]()    // c.f. format_fasta_benchmark
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

        for (size_t i = 0; i < n_entries_in_file; ++i)
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
    sequence_file_format_fastq format;

    // create temporary file, automatically removed on destruction
    test::tmp_filename file_name{"tmp.fastq"};
    auto tmp_path = file_name.get_path();

    std::ofstream ostream{tmp_path};
    ostream << generate_fastq_file(state.range(0));
    ostream.close();

    // benchmark
    for (auto _ : state)
    {
        sequence_file_input fin{tmp_path};

        // read all records and store in internal buffer
        auto it = fin.begin();
        while (it != fin.end())
            ++it;
    }
}

// ============================================================================
// seqan2 comparison
// ============================================================================
#if SEQAN3_HAS_SEQAN2

void fastq_read_from_disk_seqan2(benchmark::State & state)
{
    using namespace seqan;

    // temporary file
    test::tmp_filename file_name{"tmp.fastq"};
    auto tmp_path = file_name.get_path();

    std::ofstream ostream{tmp_path};
    ostream << generate_fastq_file(state.range(0));
    ostream.close();

    // benchmark
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
