#include <seqan3/search/algorithm/configuration/all.hpp>

int main()
{
    using namespace seqan3;
//![default]
// No errors, all hits as text position
configuration const default_cfg = search_cfg::max_error{search_cfg::total{0},
                                                        search_cfg::substitution{0},
                                                        search_cfg::insertion{0},
                                                        search_cfg::deletion{0}}
                                  | search_cfg::output{search_cfg::text_position}
                                  | search_cfg::mode{search_cfg::all};
//![default]
    return 0;
}
