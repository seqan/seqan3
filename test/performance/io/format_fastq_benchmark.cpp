// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <cstring>
#include <seqan3/std/filesystem>
#include <fstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/test/tmp_filename.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/seq_io.h>
#endif // SEQAN3_HAS_SEQAN2

constexpr unsigned default_seed = 1234u;

// ============================================================================
// generate fastq file
// ============================================================================

constexpr size_t default_sequence_length = 50; //length of nucleotide and quality sequence

std::string generate_fastq_file(size_t const entries_size)
{
    std::ostringstream stream_buffer{};
    seqan3::debug_stream_type format_stream{stream_buffer};

    std::string const id{"@name"};

    auto seed = default_seed;

    for (size_t i = 0; i < entries_size; ++i, ++seed)
    {
        auto random_sequence = seqan3::test::generate_sequence<seqan3::dna5>(default_sequence_length, 0, seed);
        auto random_qualities = seqan3::test::generate_sequence<seqan3::phred42>(default_sequence_length, 0, seed);
        format_stream << id << '\n'  // write id
                      << random_sequence << '\n' // write sequence
                      << '+' << '\n' // write qualities separator
                      << random_qualities << '\n'; // write qualities
    }

    return stream_buffer.str();
}

// ============================================================================
// save file on disc temporarily
// ============================================================================

auto create_fastq_file(size_t const entries_size)
{
    // create temporary file, automatically removed on destruction
    seqan3::test::tmp_filename fastq_file{"format_fastq_benchmark_test_file.fastq"};
    auto fastq_file_path = fastq_file.get_path();

    // fill temporary file with a fastq file
    std::ofstream ostream{fastq_file_path};
    ostream << generate_fastq_file(entries_size);
    ostream.close();
    return fastq_file;
}

// ============================================================================
// write 3-line entry to stream as often as possible
// ============================================================================

void fastq_write_to_stream_seqan3(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan3::sequence_file_output fout{ostream, seqan3::format_fastq{}};

    std::string const id{"@name"};
    auto seq = seqan3::test::generate_sequence<seqan3::dna5>(default_sequence_length, 0, default_seed);
    auto qual = seqan3::test::generate_sequence<seqan3::phred42>(default_sequence_length, 0, default_seed);

    for (auto _ : state)
        fout.emplace_back(seq, id, qual);
}

// ============================================================================
// seqan3 benchmark
// ============================================================================

// ----------------------------------------------------------------------------
// read dummy fastq file from a stream
// ----------------------------------------------------------------------------

void fastq_read_from_stream_seqan3(benchmark::State & state)
{
    std::string file = generate_fastq_file(state.range(0));

    seqan3::dna5_vector sequence{};
    std::string id{};
    std::vector<seqan3::phred42> quality{};

    for (auto _ : state)
    {
        state.PauseTiming();
        std::istringstream istream{file};
        state.ResumeTiming();

        seqan3::sequence_file_input fin{istream, seqan3::format_fastq{}};
        for (auto && [read_sequence, read_id, read_quality] : fin)
        {
            sequence = std::move(read_sequence);
            id = std::move(read_id);
            quality = std::move(read_quality);
        }

        sequence.clear();
        id.clear();
        quality.clear();
    }
}

// ----------------------------------------------------------------------------
// read dummy fastq file from temporary file on disk
// ----------------------------------------------------------------------------

void fastq_read_from_disk_seqan3(benchmark::State & state)
{
    auto file = create_fastq_file(state.range(0));
    auto path = file.get_path();

    std::vector<seqan3::dna5_vector> sequences{};
    std::vector<std::string> ids{};
    std::vector<std::vector<seqan3::phred42>> qualities{};

    for (auto _ : state)
    {
        seqan3::sequence_file_input fin{path};

        for (auto && [read_sequence, read_id, read_quality] : fin)
        {
            sequences.push_back(std::move(read_sequence));
            ids.push_back(std::move(read_id));
            qualities.push_back(std::move(read_quality));
        }

        sequences.clear();
        ids.clear();
        qualities.clear();
    }
}

// ============================================================================
// seqan2 benchmark
// ============================================================================

#if SEQAN3_HAS_SEQAN2

// ----------------------------------------------------------------------------
// read dummy fastq file from a stream
// ----------------------------------------------------------------------------

void fastq_read_from_stream_seqan2(benchmark::State & state)
{
    std::istringstream istream{};
    size_t entries_size = state.range(0);
    auto file = generate_fastq_file(entries_size);

    auto restart_iterator = [&istream, &file]()    // cf. format_fasta_benchmark
    {
        istream = std::istringstream{file};
        // same constant as seqan3 benchmark

        seqan::VirtualStream<char, seqan::Input> comp{};
        seqan::open(comp, istream);
        return seqan::directionIterator(comp, seqan::Input());
    };

    seqan::String<char> id{};
    seqan::String<seqan::Dna5> seq{};
    seqan::String<char> qual{};

    for (auto _ : state)
    {
        state.PauseTiming();
        auto it = restart_iterator();
        state.ResumeTiming();

        for (size_t i = 0; i < entries_size; ++i)
            seqan::readRecord(id, seq, qual, it, seqan::Fastq{});

        seqan::clear(id);
        seqan::clear(seq);
        seqan::clear(qual);
    }
}

// ----------------------------------------------------------------------------
// read dummy fastq file from temporary file on disk
// ----------------------------------------------------------------------------

void fastq_read_from_disk_seqan2(benchmark::State & state)
{
    seqan3::test::tmp_filename file_name{"tmp.fastq"};
    auto tmp_path = file_name.get_path();

    std::ofstream ostream{tmp_path};
    ostream << generate_fastq_file(state.range(0));
    ostream.close();

    seqan::StringSet<seqan::String<char>> ids{};
    seqan::StringSet<seqan::String<seqan::Dna5>> seqs{};
    seqan::StringSet<seqan::String<char>> quals{};

    for (auto _ : state)
    {
        seqan::SeqFileIn seqFileIn(tmp_path.c_str());
        seqan::readRecords(ids, seqs, quals, seqFileIn);

        seqan::clear(ids);
        seqan::clear(seqs);
        seqan::clear(quals);
    }
}
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK(fastq_write_to_stream_seqan3);

BENCHMARK(fastq_read_from_stream_seqan3)->Arg(100)->Arg(1000)->Arg(10000);
BENCHMARK(fastq_read_from_disk_seqan3)->Arg(100)->Arg(1000)->Arg(10000);

#if SEQAN3_HAS_SEQAN2
BENCHMARK(fastq_read_from_stream_seqan2)->Arg(100)->Arg(1000)->Arg(10000);
BENCHMARK(fastq_read_from_disk_seqan2)->Arg(100)->Arg(1000)->Arg(10000);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK_MAIN();
