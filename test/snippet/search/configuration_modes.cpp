#include <seqan3/search/algorithm/configuration/all.hpp>

int main()
{
    // Report the hit with the least number of errors (either 0 or 1 errors).
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                     seqan3::search_cfg::substitution{0},
                                                                     seqan3::search_cfg::insertion{1},
                                                                     seqan3::search_cfg::deletion{1}} |
                                       seqan3::search_cfg::mode{seqan3::search_cfg::best};

    // Report all hits with best + 1 error. E.g., if the best hit has 1 error, all hits with 1 and 2 error are reported.
    seqan3::configuration const cfg2 = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{0},
                                                                     seqan3::search_cfg::insertion{1},
                                                                     seqan3::search_cfg::deletion{1}} |
                                       seqan3::search_cfg::mode{seqan3::search_cfg::strata{1}};
    return 0;
}
