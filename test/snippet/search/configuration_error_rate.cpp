#include <seqan3/search/algorithm/configuration/max_error_rate.hpp>

int main()
{
    using namespace seqan3;
    // Allow 10% errors of any type.
    configuration const cfg1 = search_cfg::max_error_rate{search_cfg::total{0.1}};

    // Do not allow substitutions. Allow at most 10% errors.
    configuration const cfg2 = search_cfg::max_error_rate{search_cfg::total{0.1},
                                                          search_cfg::substitution{0},
                                                          search_cfg::insertion{0.1},
                                                          search_cfg::deletion{0.1}};

    // Sets total errors to 20%.
    configuration const cfg3 = search_cfg::max_error_rate{search_cfg::substitution{0},
                                                          search_cfg::insertion{0.1},
                                                          search_cfg::deletion{0.1}};

    return 0;
}
