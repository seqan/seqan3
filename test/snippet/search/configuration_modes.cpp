#include <seqan3/search/algorithm/configuration/all.hpp>

int main()
{
    using namespace seqan3;

    // Report the hit with the least number of errors (either 0 or 1 errors).
    configuration const cfg1 = search_cfg::max_error{search_cfg::total{1},
                                                     search_cfg::substitution{0},
                                                     search_cfg::insertion{1},
                                                     search_cfg::deletion{1}}
                               | search_cfg::mode{search_cfg::best};

    // Report all hits with best + 1 error. E.g., if the best hit has 1 error, all hits with 1 and 2 error are reported.
    configuration const cfg2 = search_cfg::max_error{search_cfg::substitution{0},
                                                     search_cfg::insertion{1},
                                                     search_cfg::deletion{1}}
                               | search_cfg::mode{search_cfg::strata{1}};
    return 0;
}
