#include <seqan3/search/algorithm/configuration/all.hpp>

int main()
{
    // Report hits as positions in the text.
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                     seqan3::search_cfg::substitution{0},
                                                                     seqan3::search_cfg::insertion{1},
                                                                     seqan3::search_cfg::deletion{1}} |
                                       seqan3::search_cfg::output{seqan3::search_cfg::text_position};

    // Return cursors of the index.
    seqan3::configuration const cfg2 = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{0},
                                                                     seqan3::search_cfg::insertion{1},
                                                                     seqan3::search_cfg::deletion{1}} |
                                       seqan3::search_cfg::output{seqan3::search_cfg::index_cursor};
    return 0;
}
