// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#if __has_include(<seqan/seq_io.h>)
    #include <seqan/seq_io.h>
#endif

#include <algorithm>
#include <cctype>
#include <cstring>

#include <benchmark/benchmark.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/test/performance/units.hpp>

#include <sstream>

inline constexpr size_t iterations_per_run = 1024;

inline std::string const fasta_hdr{"seq foobar blobber"};
inline std::string const fasta_seq{
    "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
    "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
    "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
    "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
    "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
    "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"};
//TODO: benchmark with spaces/newlines

static std::string fasta_file = []()
{
    std::string file{};
    for (size_t idx = 0; idx < iterations_per_run; idx++)
        file += "> " + fasta_hdr + "\n" + fasta_seq + "\n";
    return file;
}();

void write3(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan3::sequence_file_output fout{ostream, seqan3::format_fasta{}, seqan3::fields<seqan3::field::id,
                                                                                      seqan3::field::seq>{}};

    for (auto _ : state)
    {
        for (size_t i = 0; i < iterations_per_run; ++i)
            fout.emplace_back(fasta_seq, fasta_hdr);
    }

    ostream = std::ostringstream{};
    fout.emplace_back(fasta_seq, fasta_hdr);
    size_t bytes_per_run = ostream.str().size() * iterations_per_run;
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

BENCHMARK(write3);

#if __has_include(<seqan/seq_io.h>)

void write2(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan::CharString id = fasta_hdr;
    seqan::Dna5String seq = fasta_seq;

    for (auto _ : state)
    {
        for (size_t i = 0; i < iterations_per_run; ++i)
            seqan::writeRecord(ostream, id, seq, seqan::Fasta());
    }

    ostream = std::ostringstream{};
    seqan::writeRecord(ostream, id, seq, seqan::Fasta());
    size_t bytes_per_run = ostream.str().size() * iterations_per_run;
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

BENCHMARK(write2);
#endif

void read3(benchmark::State & state)
{
    std::istringstream istream{fasta_file};
    seqan3::sequence_file_input fin{istream, seqan3::format_fasta{}};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);

        auto it = fin.begin();
        for (size_t i = 0; i < iterations_per_run; ++i)
            it++;
    }

    size_t bytes_per_run = fasta_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(read3);

#if __has_include(<seqan/seq_io.h>)

#include <fstream>

void read2(benchmark::State & state)
{
    seqan::CharString id;
    seqan::Dna5String seq;

    std::istringstream istream{fasta_file};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);
        auto it = seqan::Iter<std::istringstream, seqan::StreamIterator<seqan::Input> >(istream);

        for (size_t i = 0; i < iterations_per_run; ++i)
        {
            readRecord(id, seq, it, seqan::Fasta{});
            clear(id);
            clear(seq);
        }
    }

    size_t bytes_per_run = fasta_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(read2);
#endif

BENCHMARK_MAIN();
