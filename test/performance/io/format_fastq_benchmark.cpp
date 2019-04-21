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

using namespace seqan3;

std::size_t POWER = 5; //files will have |10^POWER| lines

auto const DNA_SEQ = "AGCTAGCAGCGATCGCGATCGATCAGCGATCGAGGAATATAT"_dna5;
auto const QUALITY = "IIIIIHIIIIIIIIIIIIIIIIIIIIIIHHGIIIIHHGIIIH"_phred42;

std::size_t const iterations_per_run = 4096;  // arbitrary number directly taken from fasta benchmark

// ============================================================================
// generate dummy fastq file with 10^power identical 4-line entries
// ============================================================================
std::string generate_dummy_fastq_file(std::size_t power)
{
    std::string file{};

    std::string const id{"@name"};
    std::string const seq = DNA_SEQ | view::to_char;
    std::string const quality = QUALITY | view::to_char;

    auto const n_iterations = static_cast<unsigned int>(pow(10, power));

    for (size_t i = 0; i < n_iterations; ++i)
    {
        file += id + '\n' + seq + '\n' + '+' + '\n' + quality + '\n';
    }
    
    return file;
}

// ============================================================================
// try to write 4-line entry to stream as often as possible
// ============================================================================
void fastq_write(benchmark::State & state)
{
    std::ostringstream ostream;

    sequence_file_format_fastq format;
    sequence_file_output_options options{};

    std::string const id{"@name"};

    for (auto _ : state)
    {
        for (size_t i = 0; i < iterations_per_run; ++i)
            format.write(ostream, options, DNA_SEQ, id, QUALITY);
    }
}

BENCHMARK(fastq_write);

// ============================================================================
// read dummy fastq file with all parameters
// ============================================================================
void fastq_read_with_quality(benchmark::State & state)
{
    auto const fastq_file = generate_dummy_fastq_file(POWER);

    sequence_file_format_fastq format;
    sequence_file_input_options<dna5, false> const options{};

    std::string id{};
    std::vector<dna5> seq{};
    std::vector<phred42> qual{};

    for (auto _ : state)
    {
        state.PauseTiming();
        std::istringstream istream{fastq_file}; // costly so pause/resume is necessary
        state.ResumeTiming();

        for (size_t i = 0; i < iterations_per_run; ++i)
        {
            format.read(istream, options, seq, id, qual);
            id.clear();
            seq.clear();
            qual.clear();
        }
    }
}

BENCHMARK(fastq_read_with_quality);

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

BENCHMARK_MAIN();
