#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/parallel.hpp>

int main()
{
    // Enable parallel execution of the search algorithm with 8 threads (and allow 1 error of any type).
    seqan3::configuration const cfg1 = seqan3::search_cfg::parallel{8} |
                                       seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}};

    return 0;
}
