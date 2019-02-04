// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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
#include <seqan3/range/view/convert.hpp>

#include <sstream>

using namespace seqan3;

void write3(benchmark::State & state)
{
    std::ostringstream ostream;
    sequence_file_format_fasta format;
    sequence_file_output_options options;
    std::string id{"seq"};
    dna5_vector seq{"ACTAGACTAGCTACGATCAGCTACGATCAGCTACGA"_dna5};

    for (auto _ : state)
    {
        format.write(ostream, options, seq, id, std::ignore);
    }
}

BENCHMARK(write3);

#if __has_include(<seqan/seq_io.h>)

void write2(benchmark::State & state)
{
    std::ostringstream outStream;
    seqan::CharString meta = "seq";
    seqan::Dna5String seq = "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGA";

    for (auto _ : state)
    {
        seqan::writeRecord(outStream, meta, seq, seqan::Fasta());
    }
}

BENCHMARK(write2);
#endif

void read3(benchmark::State & state)
{
    std::string dummy_file{};
    for (size_t idx = 0; idx < 10000000; idx++)
        dummy_file += ">seq\nACTAGACTAGCTACGATCAGCTACGATCAGCTACGA\n";
    std::istringstream istream{dummy_file};
    sequence_file_format_fasta format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    dna5_vector seq;
    for (auto _ : state)
    {
        format.read(istream, options, seq, id, std::ignore);
	    id.clear();
	    seq.clear();
    }
}
BENCHMARK(read3);

#if __has_include(<seqan/seq_io.h>)

#include <fstream>

void read2(benchmark::State & state)
{
    seqan::CharString meta;
    seqan::Dna5String seq;
    std::string dummy_file{};
    for (size_t idx = 0; idx < 10000000; idx++)
        dummy_file += ">seq\nACTAGACTAGCTACGATCAGCTACGATCAGCTACGA\n";
    std::istringstream istream{dummy_file};

    seqan::VirtualStream<char, seqan::Input> comp;
    open(comp, istream);
    auto it = seqan::directionIterator(comp, seqan::Input());

    for (auto _ : state)
    {
        readRecord(meta, seq,  it, seqan::Fasta{});
        clear(meta);
        clear(seq);
    }
}
BENCHMARK(read2);
#endif

BENCHMARK_MAIN();
