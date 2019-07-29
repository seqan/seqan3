#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/std/filesystem>

int main()
{
    using seqan3::operator""_dna4;

    auto tmp_file = std::filesystem::temp_directory_path() / "my.fasta";

    {
        // Create a /tmp/my.fasta file.
        seqan3::sequence_file_output fout{tmp_file};
        fout.emplace_back("ACGT"_dna4, "TEST1");
        fout.emplace_back("AGGCTGA"_dna4, "Test2");
        fout.emplace_back("GGAGTATAATATATATATATATAT"_dna4, "Test3");
    }

    // FastA with DNA sequences assumed, regular std::ifstream taken as stream
    seqan3::sequence_file_input fin{tmp_file};

    std::filesystem::remove(tmp_file);
}
