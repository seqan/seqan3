// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/tmp_directory.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/bam_io.h>
#endif

// ============================================================================
// generate sam file from randomly generated sequence and store in file_dict
// ============================================================================

std::map<std::size_t, std::string> file_dict;

static std::string create_sam_file_string(size_t const n_queries)
{
    if (file_dict.find(n_queries) == file_dict.end())
    {
        size_t const seed{1234u};
        size_t const length_variance{0u};
        size_t const reference_size{500u};
        size_t const read_size{100u}; // typical illumina read
        std::string const query_prefix{"query_"};
        std::string const reference_id{"reference_id"};

        // generate sequences
        auto reference = seqan3::test::generate_sequence<seqan3::dna4>(reference_size, length_variance, seed);

        // align
        constexpr auto nt_score_scheme =
            seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-2}};
        auto config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                       seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                       seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                       seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
                    | seqan3::align_cfg::scoring_scheme{nt_score_scheme}
                    | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                         seqan3::align_cfg::extension_score{-1}}
                    | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_alignment{}
                    | seqan3::align_cfg::output_score{};

        using sam_fields = seqan3::fields<seqan3::field::seq,
                                          seqan3::field::id,
                                          seqan3::field::ref_id,
                                          seqan3::field::ref_offset,
                                          seqan3::field::cigar,
                                          seqan3::field::mapq,
                                          seqan3::field::qual,
                                          seqan3::field::flag>;
        std::ostringstream stream;
        seqan3::sam_file_output sam_out{stream, seqan3::format_sam{}, sam_fields{}};

        for (size_t i = 0; i < n_queries; ++i)
        {
            auto query = seqan3::test::generate_sequence<seqan3::dna4>(read_size, length_variance, seed + i);
            auto qualities = seqan3::test::generate_sequence<seqan3::phred42>(read_size, length_variance, seed + i);
            auto align_result = *(seqan3::align_pairwise(std::tie(reference, query), config).begin());
            std::string const current_query_id = query_prefix + std::to_string(i);

            sam_out.emplace_back(query,                                                  // field::seq
                                 current_query_id,                                       // field::id
                                 reference_id,                                           // field::ref_id
                                 align_result.sequence1_begin_position(),                // field::ref_offset
                                 seqan3::cigar_from_alignment(align_result.alignment()), // field::cigar
                                 align_result.score(),                                   // field::mapq
                                 qualities,                                              // field::qual
                                 seqan3::sam_flag::none);                                // field::flag
        }

        file_dict[n_queries] = stream.str();
    }

    return file_dict[n_queries];
}

void write_file(std::string const & file_name, size_t const n_queries)
{
    std::ofstream ostream{file_name};
    ostream << create_sam_file_string(n_queries);
    ostream.close();
}

// ============================================================================
// seqan3
// ============================================================================

void sam_file_read_from_stream(benchmark::State & state)
{
    size_t const n_queries = state.range(0);

    std::istringstream istream{create_sam_file_string(n_queries)};

    for (auto _ : state)
    {
        // refill stream (same constant for seqan2 benchmark)
        istream.clear();
        istream.seekg(0, std::ios::beg);

        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};

        // read all records and store in internal buffer
        auto it = fin.begin();
        while (it != fin.end())
            ++it;
    }
}

void sam_file_read_from_disk(benchmark::State & state)
{
    size_t const n_queries = state.range(0);

    seqan3::test::tmp_directory tmp{};
    auto tmp_path = tmp.path() / "tmp.sam";

    write_file(tmp_path, n_queries);

    for (auto _ : state)
    {
        seqan3::sam_file_input fin{tmp_path};

        // read all records and store in internal buffer
        auto it = fin.begin();
        while (it != fin.end())
            ++it;
    }
}

#if SEQAN3_HAS_SEQAN2
// ============================================================================
// seqan2 read from stream
// ============================================================================

void seqan2_sam_file_read_from_stream(benchmark::State & state)
{
    size_t const n_queries = state.range(0);
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "tmp.sam";

    std::string sam_file = create_sam_file_string(n_queries);

    // create temporary BamFileIn and read from disk to get the context...
    write_file(filename, n_queries);
    seqan2::BamHeader tmp_header;
    seqan2::BamFileIn tmp_bam_file_in(filename.c_str());
    seqan2::readHeader(tmp_header, tmp_bam_file_in);
    auto cxt = seqan2::context(tmp_bam_file_in);

    seqan2::BamAlignmentRecord record;
    seqan2::BamHeader header;

    std::istringstream istream{sam_file};

    for (auto _ : state)
    {
        istream.clear();
        istream.seekg(0, std::ios::beg);

        auto it = seqan2::Iter<std::istringstream, seqan2::StreamIterator<seqan2::Input>>(istream);

        seqan2::readHeader(header, cxt, it, seqan2::Sam());

        for (size_t i = 0; i < n_queries; ++i)
        {
            seqan2::readRecord(record, cxt, it, seqan2::Sam());
            seqan2::clear(record);
        }

        clear(header);
    }
}

void seqan2_sam_file_read_from_disk(benchmark::State & state)
{
    size_t const n_queries = state.range(0);

    seqan3::test::tmp_directory tmp{};
    auto tmp_path = tmp.path() / "tmp.sam";

    write_file(tmp_path, n_queries);

    seqan2::BamHeader header;
    seqan2::BamAlignmentRecord record;

    for (auto _ : state)
    {
        seqan2::BamFileIn bamFileIn(tmp_path.c_str());

        seqan2::readHeader(header, bamFileIn);

        while (!seqan2::atEnd(bamFileIn))
            seqan2::readRecord(record, bamFileIn);

        seqan2::clear(header);
        seqan2::clear(record);
    }
}

#endif // SEQAN3_HAS_SEQAN2

#ifndef NDEBUG
static constexpr size_t low_query_count{5u};
static constexpr size_t high_query_count{10u};
#else
static constexpr size_t low_query_count{50u};
static constexpr size_t high_query_count{500u};
#endif // NDEBUG

BENCHMARK(sam_file_read_from_stream)->Arg(low_query_count);
BENCHMARK(sam_file_read_from_stream)->Arg(high_query_count);

BENCHMARK(sam_file_read_from_disk)->Arg(low_query_count);
BENCHMARK(sam_file_read_from_disk)->Arg(high_query_count);

#if SEQAN3_HAS_SEQAN2
BENCHMARK(seqan2_sam_file_read_from_stream)->Arg(low_query_count);
BENCHMARK(seqan2_sam_file_read_from_stream)->Arg(high_query_count);

BENCHMARK(seqan2_sam_file_read_from_disk)->Arg(low_query_count);
BENCHMARK(seqan2_sam_file_read_from_disk)->Arg(high_query_count);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK_MAIN();
