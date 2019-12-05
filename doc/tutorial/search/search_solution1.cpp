#if SEQAN3_WITH_CEREAL

#include "cleanup.hpp"
seqan3::cleanup index_file{"index.file"};

//![solution]
#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/search/fm_index/fm_index.hpp>

using seqan3::operator""_dna4;

int main()
{
    seqan3::dna4_vector
                text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::fm_index index{text};

    {
        std::ofstream os{"index.file", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

    // we need to tell the index that we work on a single text and a `dna4` alphabet before loading
    seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single> index2;
    {
        std::ifstream is{"index.file", std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index2);
    }

    if (index == index2)
        std::cout << "The indices are identical!\n";
    else
        std::cout << "The indices differ!\n";
}
//![solution]
#endif //SEQAN3_WITH_CEREAL
