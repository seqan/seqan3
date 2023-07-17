// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <sstream>
#include <string>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/utility/range/to.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/rna_io.h>
#endif

inline constexpr size_t iterations_per_run = 1024;

inline std::string const header{"seq foobar blobber"};
inline auto const rna_sequence = seqan3::test::generate_sequence<seqan3::rna4>(474, 0, 0);
auto const sequence = rna_sequence | seqan3::views::to_char | seqan3::ranges::to<std::string>();

inline std::string const structure{"(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))......."
                                   "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))......."
                                   "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))......."
                                   "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))......."
                                   "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))......."
                                   "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))......."};

inline std::string const vienna_file = []()
{
    std::string file{};
    for (size_t idx = 0; idx < iterations_per_run; idx++)
        file += "> " + header + "\n" + sequence + "\n" + structure + "\n";
    return file;
}();

void write_seqan3(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan3::structure_file_output fout{ostream, seqan3::format_vienna{}};

    for (auto _ : state)
    {
        for (size_t idx = 0; idx < iterations_per_run; ++idx)
            fout.emplace_back(sequence, header, structure);
    }

    ostream = std::ostringstream{};
    fout.emplace_back(sequence, header, structure);
    size_t bytes_per_run = ostream.str().size() * iterations_per_run;
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(write_seqan3);

#if SEQAN3_HAS_SEQAN2
void write_seqan2(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan2::RnaRecord record{};
    record.name = header;
    record.sequence = sequence;
    seqan2::bracket2graph(record.fixedGraphs, structure);

    for (auto _ : state)
    {
        for (size_t idx = 0; idx < iterations_per_run; ++idx)
            seqan2::writeRecord(ostream, record, seqan2::Vienna());
    }

    ostream = std::ostringstream{};
    seqan2::writeRecord(ostream, record, seqan2::Vienna());
    size_t bytes_per_run = ostream.str().size() * iterations_per_run;
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(write_seqan2);
#endif

void read_seqan3(benchmark::State & state)
{
    std::istringstream istream{vienna_file};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);

        auto it = fin.begin();
        for (size_t idx = 0; idx < iterations_per_run; ++idx)
            it++;
    }

    size_t bytes_per_run = vienna_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(read_seqan3);

#if SEQAN3_HAS_SEQAN2
void read_seqan2(benchmark::State & state)
{
    seqan2::RnaRecord record{};
    std::istringstream istream{vienna_file};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);
        auto it = seqan2::Iter<std::istringstream, seqan2::StreamIterator<seqan2::Input>>(istream);

        for (size_t idx = 0; idx < iterations_per_run; ++idx)
        {
            seqan2::readRecord(record, it, seqan2::Vienna());
            clear(record);
        }
    }

    size_t bytes_per_run = vienna_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(read_seqan2);
#endif

BENCHMARK_MAIN();
