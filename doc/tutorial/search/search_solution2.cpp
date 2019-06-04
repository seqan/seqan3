#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/algorithm/search.hpp>

using namespace seqan3;

void run_text_single()
{
    dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    fm_index index{text};

    auto results = search("GCT"_dna4, index);
    std::ranges::sort(results);
    debug_stream << "=====   Running on a single text   =====\n";
    debug_stream << "There are " << results.size() << " hits.\n";
    debug_stream << "The positions are " << results << '\n';
}

void run_text_collection()
{
    std::vector<dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                  "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                  "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    fm_index index{text};

    auto results = search("GCT"_dna4, index);
    std::ranges::sort(results);
    debug_stream << "===== Running on a text collection =====\n";
    debug_stream << "There are " << results.size() << " hits.\n";
    debug_stream << "The positions are " << results << '\n';
}

int main()
{
   run_text_single();
   debug_stream << '\n';
   run_text_collection();
}
