#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sam_file/all.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<std::string> ids = {"read1", "read2"};
    std::vector<std::vector<seqan3::dna4>> seqs = {"ACGATCGACTAGCTACGATCAGCTAGCAG"_dna4,
                                                   "AGAAAGAGCGAGGCTATTTTAGCGAGTTA"_dna4};

    auto tmp_dir = std::filesystem::temp_directory_path();
    seqan3::sam_file_output fout{tmp_dir/"my.sam", seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

    for (size_t i = 0; i < ids.size(); ++i)
    {
        fout.emplace_back(ids[i], seqs[i]);
    }
}
