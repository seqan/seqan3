#include <seqan3/search/fm_index/bi_fm_index.hpp>   // for using the bi_fm_index
#include <seqan3/search/fm_index/fm_index.hpp>      // for using the fm_index

using namespace seqan3;

int main()
{
    std::string text{"Garfield the fat cat without a hat."};
    fm_index index{text};        // unidirectional index on single text
    bi_fm_index bi_index{text};  // bidirectional index on single text
}
