// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <cstring>
#include <filesystem>
#include <fstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/test/tmp_directory.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/seq_io.h>
#endif // SEQAN3_HAS_SEQAN2

constexpr unsigned default_seed = 1234u;

// ============================================================================
// generate FASTQ file
// ============================================================================

inline constexpr size_t default_sequence_length = 50; //length of nucleotide and quality sequence
inline std::string const fastq_id{"the FASTQ file"};

std::string generate_fastq_string(size_t const entries_size)
{
    std::ostringstream stream_buffer{};
    seqan3::sequence_file_output fastq_ostream{stream_buffer, seqan3::format_fastq{}};

    auto seed = default_seed;

    for (size_t i = 0; i < entries_size; ++i, ++seed)
    {
        auto random_sequence = seqan3::test::generate_sequence<seqan3::dna5>(default_sequence_length, 0, seed);
        auto random_qualities = seqan3::test::generate_sequence<seqan3::phred42>(default_sequence_length, 0, seed);

        fastq_ostream.emplace_back(random_sequence, fastq_id, random_qualities);
    }

    stream_buffer.flush();
    return stream_buffer.str();
}

// ============================================================================
// save file on disc temporarily
// ============================================================================

void create_fastq_file_for(std::filesystem::path const & fastq_file, std::string const & fastq_content)
{
    // fill temporary file with a FASTQ file
    std::ofstream ostream{fastq_file};
    ostream << fastq_content;
    ostream.close();
}

// ============================================================================
// seqan3 FASTQ output benchmark
// ============================================================================

// ----------------------------------------------------------------------------
// write dummy FASTQ file to a string stream
// ----------------------------------------------------------------------------

void fastq_write_to_stream_seqan3(benchmark::State & state)
{
    size_t const iterations_per_run = state.range(0);
    std::ostringstream ostream;
    seqan3::sequence_file_output fout{ostream,
                                      seqan3::format_fastq{},
                                      seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>{}};

    auto seq = seqan3::test::generate_sequence<seqan3::dna5>(default_sequence_length, 0, default_seed);
    auto qual = seqan3::test::generate_sequence<seqan3::phred42>(default_sequence_length, 0, default_seed);

    for (auto _ : state)
    {
        for (size_t i = 0; i < iterations_per_run; ++i)
            fout.emplace_back(fastq_id, seq, qual);
    }

    ostream.str("");                        // Reset stream.
    fout.emplace_back(fastq_id, seq, qual); // Write one entry.
    ostream.flush();                        // Make sure the buffer is flushed.

    size_t bytes_per_run = ostream.str().size() * iterations_per_run;
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

BENCHMARK(fastq_write_to_stream_seqan3)->Arg(100)->Arg(1000)->Arg(10000);

// ============================================================================
// seqan3 FASTQ input benchmark
// ============================================================================

// ----------------------------------------------------------------------------
// read dummy FASTQ file from a stream
// ----------------------------------------------------------------------------

void fastq_read_from_stream_seqan3(benchmark::State & state)
{
    size_t const iterations_per_run = state.range(0);
    std::string fastq_file = generate_fastq_string(iterations_per_run);
    std::istringstream istream{fastq_file};
    seqan3::sequence_file_input fastq_file_in{istream, seqan3::format_fastq{}};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);

        auto it = fastq_file_in.begin();
        for (size_t i = 0; i < iterations_per_run; ++i)
            it++;
    }

    size_t bytes_per_run = fastq_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

// ----------------------------------------------------------------------------
// read dummy FASTQ file from temporary file on disk
// ----------------------------------------------------------------------------

void fastq_read_from_disk_seqan3(benchmark::State & state)
{
    size_t const iterations_per_run = state.range(0);
    std::string fastq_file = generate_fastq_string(iterations_per_run);
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "format_fastq_benchmark_test_file.fastq";
    create_fastq_file_for(filename, fastq_file);

    for (auto _ : state)
    {
        seqan3::sequence_file_input fastq_file_in{filename};
        auto it = fastq_file_in.begin();
        for (size_t i = 0; i < iterations_per_run; ++i)
            it++;
    }

    size_t bytes_per_run = fastq_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

// ============================================================================
// seqan2 FASTQ input benchmark
// ============================================================================

#if SEQAN3_HAS_SEQAN2

// ----------------------------------------------------------------------------
// read dummy FASTQ file from a stream
// ----------------------------------------------------------------------------

void fastq_read_from_stream_seqan2(benchmark::State & state)
{
    size_t const iterations_per_run = state.range(0);
    std::string fastq_file = generate_fastq_string(iterations_per_run);
    std::istringstream istream{fastq_file};

    seqan2::CharString id{};
    seqan2::Dna5String seq{};
    seqan2::CharString qual{};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);
        auto it = seqan2::Iter<std::istringstream, seqan2::StreamIterator<seqan2::Input>>(istream);

        for (size_t i = 0; i < iterations_per_run; ++i)
        {
            readRecord(id, seq, qual, it, seqan2::Fastq{});
            clear(id);
            clear(seq);
            clear(qual);
        }
    }

    size_t bytes_per_run = fastq_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

// ----------------------------------------------------------------------------
// read dummy FASTQ file from temporary file on disk
// ----------------------------------------------------------------------------

void fastq_read_from_disk_seqan2(benchmark::State & state)
{
    size_t const iterations_per_run = state.range(0);
    std::string fastq_content = generate_fastq_string(iterations_per_run);
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "format_fastq_benchmark_test_file.fastq";
    create_fastq_file_for(filename, fastq_content);

    seqan2::CharString id{};
    seqan2::Dna5String seq{};
    seqan2::CharString qual{};

    for (auto _ : state)
    {
        seqan2::SeqFileIn seqFileIn(filename.c_str());
        auto it = seqFileIn.iter;

        for (size_t i = 0; i < iterations_per_run; ++i)
        {
            readRecord(id, seq, qual, it, seqan2::Fastq{});
            clear(id);
            clear(seq);
            clear(qual);
        }
    }

    size_t bytes_per_run = fastq_content.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK(fastq_read_from_stream_seqan3)->Arg(100)->Arg(1000)->Arg(10000);
BENCHMARK(fastq_read_from_disk_seqan3)->Arg(100)->Arg(1000)->Arg(10000);

#if SEQAN3_HAS_SEQAN2
BENCHMARK(fastq_read_from_stream_seqan2)->Arg(100)->Arg(1000)->Arg(10000);
BENCHMARK(fastq_read_from_disk_seqan2)->Arg(100)->Arg(1000)->Arg(10000);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK_MAIN();
