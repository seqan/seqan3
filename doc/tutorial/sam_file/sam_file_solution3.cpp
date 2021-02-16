#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/std/filesystem>

using seqan3::operator""_dna4;

int main()
{
    std::vector<std::string> ids = {"read1", "read2"};
    std::vector<std::vector<seqan3::dna4>> seqs = {"ACGATCGACTAGCTACGATCAGCTAGCAG"_dna4,
                                                   "AGAAAGAGCGAGGCTATTTTAGCGAGTTA"_dna4};

    auto tmp_dir = std::filesystem::temp_directory_path();
    seqan3::alignment_file_output fout{tmp_dir/"my.sam", seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

    for (size_t i = 0; i < ids.size(); ++i)
    {
        fout.emplace_back(ids[i], seqs[i]);
    }
}
