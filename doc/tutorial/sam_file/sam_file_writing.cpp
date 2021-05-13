#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "out.sam" will be deleted after the execution
seqan3::test::create_temporary_snippet_file example_sam{"out.sam", ""};

//![main]
#include <seqan3/io/sam_file/all.hpp>

using aligned_sequence_type = std::vector<seqan3::gapped<seqan3::dna5>>;
using alignment_type = std::pair<aligned_sequence_type, aligned_sequence_type>;

using types = seqan3::type_list<std::vector<seqan3::dna5>, std::string, alignment_type>;
using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::alignment>;

int main()
{
    using namespace seqan3::literals;

    auto filename = std::filesystem::current_path() / "out.sam";

    seqan3::sam_file_output fout{filename};
    using sam_record_type = seqan3::sam_record<types, fields>;

    // write the following to the file
    // r001	0	*	0	0	4M2I2M2D	*	0	0	ACGTACGT	*
    sam_record_type record{};
    record.id() = "r001";
    record.sequence() = "ACGTACGT"_dna5;
    auto & [reference_sequence, read_sequence] = record.alignment();

    // ACGT--GTTT
    seqan3::assign_unaligned(reference_sequence, "ACGTGTTT"_dna5);
    seqan3::insert_gap(reference_sequence, reference_sequence.begin() + 4, 2);

    // ACGTACGT--
    seqan3::assign_unaligned(read_sequence, record.sequence());
    seqan3::insert_gap(read_sequence, read_sequence.end(), 2);

    fout.push_back(record);
}
//![main]
