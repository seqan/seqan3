#include <seqan3/search/configuration/all.hpp>

int main()
{
    seqan3::search_cfg::mode my_mode{}; // Enables the dynamic configuration, defaults to seqan3::search_cfg::all.

    if (true)
        my_mode = seqan3::search_cfg::best; // The best mode is selected.
    else
        my_mode = seqan3::search_cfg::strata{2}; // The strata mode is selected.

    // Combine the dynamically configured mode with other search configurations.
    seqan3::configuration const cfg = my_mode |
                                      seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                    seqan3::search_cfg::substitution{0},
                                                                    seqan3::search_cfg::insertion{1},
                                                                    seqan3::search_cfg::deletion{1}};
    return 0;
}
