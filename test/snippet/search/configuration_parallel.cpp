#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/parallel.hpp>

int main()
{
    // Enable parallel execution of the search algorithm with 8 threads (and allow 1 error of any type).
    seqan3::configuration const cfg1 = seqan3::search_cfg::parallel{8} |
                                       seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    return 0;
}
