// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <sstream>
#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/io/sequence_file/format_embl.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/utility/range/to.hpp>

inline constexpr size_t iterations_per_run = 1024;

inline std::string const embl_id{"seq foobar blobber"};
// 494 bp: 8 full lines of 60 bp and a last line of 14 bp.
inline seqan3::dna5_vector const embl_seq =
    std::string{"ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
                "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
                "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
                "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
                "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
                "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGAACTAGACTAGCTACGATACTAGACTAGCTACGATCAGCTACGA"
                "ACTAGACTAGCTACGATCAG"}
    | seqan3::views::char_to<seqan3::dna5> | seqan3::ranges::to<seqan3::dna5_vector>();

static std::string embl_file = []()
{
    std::ostringstream ostream;

    {
        seqan3::sequence_file_output fout{ostream,
                                          seqan3::format_embl{},
                                          seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

        for (size_t idx = 0; idx < iterations_per_run; idx++)
            fout.emplace_back(embl_id, embl_seq);
    }

    return ostream.str();
}();

void seqan3_dna5_ostringstream_write(benchmark::State & state)
{
    std::ostringstream ostream;
    seqan3::sequence_file_output fout{ostream,
                                      seqan3::format_embl{},
                                      seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

    for (auto _ : state)
    {
        for (size_t i = 0; i < iterations_per_run; ++i)
            fout.emplace_back(embl_id, embl_seq);
    }

    ostream = std::ostringstream{};
    fout.emplace_back(embl_id, embl_seq);
    size_t bytes_per_run = ostream.str().size() * iterations_per_run;
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}

BENCHMARK(seqan3_dna5_ostringstream_write);

void seqan3_dna5_istringstream_read(benchmark::State & state)
{
    std::istringstream istream{embl_file};
    seqan3::sequence_file_input fin{istream, seqan3::format_embl{}};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);

        auto it = fin.begin();
        for (size_t i = 0; i < iterations_per_run; ++i)
            it++;
    }

    size_t bytes_per_run = embl_file.size();
    state.counters["iterations_per_run"] = iterations_per_run;
    state.counters["bytes_per_run"] = bytes_per_run;
    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(bytes_per_run);
}
BENCHMARK(seqan3_dna5_istringstream_read);

BENCHMARK_MAIN();
