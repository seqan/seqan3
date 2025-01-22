// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#if __has_include(<seqan/seq_io.h>)
#    include <seqan/seq_io.h>
#endif

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <sstream>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/utility/views/convert.hpp>

inline constexpr size_t iterations_per_run = 1024;

inline std::string const fasta_hdr{"seq foobar blobber"};
inline std::string const fasta_seq{"ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
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

void seqan3_dna5_ostringstream_write(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan3::sequence_file_output fout{ostream,
                                      seqan3::format_fasta{},
                                      seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

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

BENCHMARK(seqan3_dna5_ostringstream_write);

#if __has_include(<seqan/seq_io.h>)

void seqan2_dna5_ostringstream_write(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan2::CharString id = fasta_hdr;
    seqan2::Dna5String seq = fasta_seq;

    for (auto _ : state)
    {
        for (size_t i = 0; i < iterations_per_run; ++i)
            seqan2::writeRecord(ostream, id, seq, seqan2::Fasta());
    }

    ostream = std::ostringstream{};
    seqan2::writeRecord(ostream, id, seq, seqan2::Fasta());
    size_t bytes_per_run = ostream.str().size() * iterations_per_run;
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

BENCHMARK(seqan2_dna5_ostringstream_write);
#endif

void seqan3_dna5_istringstream_read(benchmark::State & state)
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
BENCHMARK(seqan3_dna5_istringstream_read);

#if __has_include(<seqan/seq_io.h>)

#    include <fstream>

void seqan2_dna5_istringstream_read(benchmark::State & state)
{
    seqan2::CharString id;
    seqan2::Dna5String seq;

    std::istringstream istream{fasta_file};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);
        auto it = seqan2::Iter<std::istringstream, seqan2::StreamIterator<seqan2::Input>>(istream);

        for (size_t i = 0; i < iterations_per_run; ++i)
        {
            readRecord(id, seq, it, seqan2::Fasta{});
            clear(id);
            clear(seq);
        }
    }

    size_t bytes_per_run = fasta_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(seqan2_dna5_istringstream_read);
#endif

BENCHMARK_MAIN();
