#include "cleanup.hpp"
seqan3::cleanup index_file{"index.file"};

#if SEQAN3_WITH_CEREAL
#include <fstream>

#include <cereal/archives/binary.hpp>
#endif //SEQAN3_WITH_CEREAL

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/filesystem>

using namespace seqan3;
using namespace std::string_literals;

int main()
{

{
//![text_collection]
std::vector<std::string> text{{"Garfield the fat cat without a hat."},
                              {"He is infinite, he is eternal."},
                              {"Yet another text I have to think about."}};
fm_index index{text};
bi_fm_index bi_index{text};
//![text_collection]
}

#if SEQAN3_WITH_CEREAL
{
//![store]
#include <fstream>                                  // for writing/reading files

#include <cereal/archives/binary.hpp>               // for storing/loading indices via cereal

std::string text{"Garfield the fat cat without a hat."};
fm_index index{text};
{
    std::ofstream os{"index.file", std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}
//![store]
}

{
//![load]
fm_index<text_layout::single> index; // we need to tell the index that we work on a single text before loading
{
    std::ifstream is{"index.file", std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(index);
}
//![load]
}
#endif //SEQAN3_WITH_CEREAL

{
//![error_search]
std::string text{"Garfield the fat cat without a hat."};
fm_index index{text};
configuration const cfg = search_cfg::max_error{search_cfg::total{1},
                                                search_cfg::substitution{0},
                                                search_cfg::insertion{1},
                                                search_cfg::deletion{1}};
debug_stream << search("cat"s, index, cfg) << '\n'; // [14,17,18,32]
//![error_search]
}

{
//![multiple_queries]
std::string text{"Garfield the fat cat without a hat."};
fm_index index{text};
std::vector<std::string> query{"cat"s, "hat"s};
debug_stream << search(query, index) << '\n'; // [[17],[31]]
//![multiple_queries]
}

{
//![error_sum]
configuration const cfg = search_cfg::max_error{search_cfg::total{2},
                                                search_cfg::substitution{2},
                                                search_cfg::insertion{1},
                                                search_cfg::deletion{1}};
//![error_sum]
}

{
//![mode_best]
configuration const cfg = search_cfg::max_error{search_cfg::total{1},
                                                search_cfg::substitution{0},
                                                search_cfg::insertion{1},
                                                search_cfg::deletion{1}} |
                          search_cfg::mode{search_cfg::best};
//![mode_best]
}

{
//![mode_strata]
configuration const cfg = search_cfg::max_error{search_cfg::total{2},
                                                search_cfg::substitution{0},
                                                search_cfg::insertion{1},
                                                search_cfg::deletion{1}} |
                          search_cfg::mode{search_cfg::strata{2}};
//![mode_strata]
}

}
