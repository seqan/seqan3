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
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

using namespace seqan3;

unsigned int const SEED = 1234;

// ============================================================================
// generate new file
// ============================================================================
static constexpr std::size_t SEQUENCE_LENGTH = 300;
static constexpr std::size_t N_ENTRIES_IN_FILE = 4069; // number of 3-line entries

// workaround because test::generate_sequence does not work at compile time
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
        return file;
    }
    else
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
// read dummy fastq file with all parameters
// ============================================================================
void fastq_read_with_quality(benchmark::State & state)
{
    sequence_file_format_fastq format;
    sequence_file_input_options<dna5, false> const options{};

    std::string id{};
    std::vector<dna5> seq{};
    std::vector<phred42> qual{};

    std::istringstream original_istream{get_file()};
    std::istringstream swap_istream{get_file};

    for (auto _ : state)
    {
        format.read(swap_istream, options, seq, id, qual);

        // all O(1)
        id.clear();
        seq.clear();
        qual.clear();

        // refill stream

    }
}

BENCHMARK(fastq_read_with_quality);

/*
// ============================================================================
// read dummy fastq file ignoring only quality
// ============================================================================
void fastq_read_without_quality(benchmark::State & state)
{
    auto const fastq_file = generate_dummy_fastq_file(POWER);

    sequence_file_format_fastq format;
    sequence_file_input_options<dna5, false> const options{};

    std::string id{};
    std::vector<dna5> seq{};

    for (auto _ : state)
    {
        state.PauseTiming();
        std::istringstream istream{fastq_file};
        state.ResumeTiming();

        for (size_t i = 0; i < iterations_per_run; ++i)
        {
            format.read(istream, options, seq, id, std::ignore);
            id.clear();
            seq.clear();
        }
    }
}

BENCHMARK(fastq_read_without_quality);

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

// TODO <Clemens C.>: seqan2 comparison
// TODO <Clemens C.>: refactor fasta benchmark for exact comparison
*/
BENCHMARK_MAIN();
