#include <sstream>
#include <string>
#include <tuple>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

using seqan3::operator""_dna4;

struct data_storage_t
{
    seqan3::concatenated_sequences<seqan3::dna4_vector> sequences{"ACGT"_dna4, "AAA"_dna4};
    seqan3::concatenated_sequences<std::string> ids{std::string{"ID1"}, std::string{"ID2"}};
};

int main()
{
    data_storage_t data_storage{};

    // ... in your file writing function:

    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};

    fout = std::tie(data_storage.sequences, data_storage.ids);
}
