#if SEQAN3_WITH_CEREAL

#include "cleanup.hpp"
seqan3::cleanup index_file{"index.file"};

//![solution]
#include <fstream>

#include <cereal/archives/binary.hpp>

#include <seqan3/search/fm_index/fm_index.hpp>

using namespace seqan3;

int main()
{
    dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    fm_index index{text};

    {
        std::ofstream os{"index.file", std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

    fm_index<text_layout::single> index2; // we need to tell the index that we work on a single text before loading
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
