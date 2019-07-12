#include <seqan3/search/algorithm/configuration/max_error_rate.hpp>

int main()
{
    // Allow 10% errors of any type.
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{0.1}};

    // Do not allow substitutions. Allow at most 10% errors.
    seqan3::configuration const cfg2 = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{0.1},
                                                                          seqan3::search_cfg::substitution{0},
                                                                          seqan3::search_cfg::insertion{0.1},
                                                                          seqan3::search_cfg::deletion{0.1}};

    // Sets total errors to 20%.
    seqan3::configuration const cfg3 = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::substitution{0},
                                                                          seqan3::search_cfg::insertion{0.1},
                                                                          seqan3::search_cfg::deletion{0.1}};

    return 0;
}
