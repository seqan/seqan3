#include <seqan3/search/algorithm/configuration/max_error.hpp>

int main()
{
    // Allow 1 error of any type.
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}};

    // Do not allow substitutions. Allow at most 1 error.
    seqan3::configuration const cfg2 = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                     seqan3::search_cfg::substitution{0},
                                                                     seqan3::search_cfg::insertion{1},
                                                                     seqan3::search_cfg::deletion{1}};

    // Sets total errors to 2.
    seqan3::configuration const cfg3 = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{0},
                                                                     seqan3::search_cfg::insertion{1},
                                                                     seqan3::search_cfg::deletion{1}};

    return 0;
}
