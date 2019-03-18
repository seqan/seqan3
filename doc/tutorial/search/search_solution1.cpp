#include "cleanup.hpp"
seqan3::cleanup index_file{"index.file"};
seqan3::cleanup index_file_tb{"index.file.tb"};
seqan3::cleanup index_file_tbrs{"index.file.tbrs"};
seqan3::cleanup index_file_tbss{"index.file.tbss"};

//![solution]
#include <seqan3/search/fm_index/fm_index.hpp>

using namespace seqan3;

int main()
{
    dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    fm_index index{text};

    if (index.store("index.file"))
        std::cout << "Index stored successfully.\n";
    else
        std::cout << "Index could not be stored.\n";

    fm_index<dna4_vector> index2;

    if (index2.load("index.file"))
        std::cout << "Index loaded successfully.\n";
    else
        std::cout << "Index could not be loaded.\n";
}
//![solution]
