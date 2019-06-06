#include <seqan3/search/algorithm/configuration/max_error.hpp>

int main()
{
    using namespace seqan3;
    // Allow 1 error of any type.
    configuration const cfg1 = search_cfg::max_error{search_cfg::total{1}};

    // Do not allow substitutions. Allow at most 1 error.
    configuration const cfg2 = search_cfg::max_error{search_cfg::total{1},
                                                     search_cfg::substitution{0},
                                                     search_cfg::insertion{1},
                                                     search_cfg::deletion{1}};

    // Sets total errors to 2.
    configuration const cfg3 = search_cfg::max_error{search_cfg::substitution{0},
                                                     search_cfg::insertion{1},
                                                     search_cfg::deletion{1}};

    return 0;
}
