#include <seqan3/search/algorithm/configuration/all.hpp>

int main()
{
    using namespace seqan3;

    // Report hits as positions in the text.
    configuration const cfg1 = search_cfg::max_error{search_cfg::total{1},
                                                     search_cfg::substitution{0},
                                                     search_cfg::insertion{1},
                                                     search_cfg::deletion{1}}
                               | search_cfg::output{search_cfg::text_position};

    // Return cursors of the index.
    configuration const cfg2 = search_cfg::max_error{search_cfg::substitution{0},
                                                     search_cfg::insertion{1},
                                                     search_cfg::deletion{1}}
                               | search_cfg::output{search_cfg::index_cursor};
    return 0;
}
