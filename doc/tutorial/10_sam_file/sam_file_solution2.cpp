// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file reference_fasta{"reference.fasta",
                                                            R"//![ref_file](
>chr1
ACAGCAGGCATCTATCGGCGGATCGATCAGGCAGGCAGCTACTGG
>chr2
ACAGCAGGCATCTATCGGCGGATCGATCAGGCAGGCAGCTACTGTAATGGCATCAAAATCGGCATG
)//![ref_file]"}; // std::filesystem::current_path() / "reference.fasta" will be deleted after the execution

seqan3::test::create_temporary_snippet_file mapping_sam{"mapping.sam",
                                                        R"//![sam_file](
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:45
@SQ	SN:chr2	LN:66
r001	99	chr1	7	60	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	chr1	9	60	5S6M	*	0	0	GCCTAAGCTAA	*
r004	0	chr2	16	60	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	chr2	18	10	5M	*	0	0	TAGGC	*
)//![sam_file]"}; // std::filesystem::current_path() / "mapping.sam" will be deleted after the execution

//![solution]
#include <algorithm> // std::ranges::count
#include <filesystem>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sequence_file/input.hpp>

int main()
{
    std::filesystem::path current_path = std::filesystem::current_path();

    // read in reference information
    seqan3::sequence_file_input reference_file{current_path / "reference.fasta"};

    std::vector<std::string> reference_ids{};
    std::vector<seqan3::dna5_vector> reference_sequences{};

    for (auto && record : reference_file)
    {
        reference_ids.push_back(std::move(record.id()));
        reference_sequences.push_back(std::move(record.sequence()));
    }

    // filter out alignments
    seqan3::sam_file_input mapping_file{current_path / "mapping.sam", reference_ids, reference_sequences};

    auto mapq_filter = std::views::filter(
        [](auto & record)
        {
            return record.mapping_quality() >= 30;
        });

    for (auto & record : mapping_file | mapq_filter)
    {
        auto alignment = seqan3::alignment_from_cigar(record.cigar_sequence(),
                                                      reference_sequences[record.reference_id().value()],
                                                      record.reference_position().value(),
                                                      record.sequence());

        // as loop
        size_t sum_reference{};
        for (auto const & char_reference : std::get<0>(alignment))
            if (char_reference == seqan3::gap{})
                ++sum_reference;

        // or via std::ranges::count
        size_t sum_read = std::ranges::count(std::get<1>(alignment), seqan3::gap{});

        // The reference_id is ZERO based and an optional. -1 is represented by std::nullopt (= reference not known).
        std::optional reference_id = record.reference_id();

        seqan3::debug_stream << record.id() << " mapped against "
                             << (reference_id ? std::to_string(reference_id.value()) : "unknown reference") << " with "
                             << sum_read << " gaps in the read sequence and " << sum_reference
                             << " gaps in the reference sequence.\n";
    }
}
//![solution]
