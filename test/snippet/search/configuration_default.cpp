#include <seqan3/search/algorithm/configuration/all.hpp>

int main()
{
    // No errors, all hits as text position
    seqan3::configuration const default_cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0},
                                                                            seqan3::search_cfg::substitution{0},
                                                                            seqan3::search_cfg::insertion{0},
                                                                            seqan3::search_cfg::deletion{0}} |
                                              seqan3::search_cfg::output{seqan3::search_cfg::text_position} | 
                                              seqan3::search_cfg::mode{seqan3::search_cfg::all};
    return 0;
}
