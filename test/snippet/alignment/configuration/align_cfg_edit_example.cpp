#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_max_error.hpp>

int main()
{
{
//! [semi_global]
    using namespace seqan3;

    // Computes semi global edit distance using fast-bit vector algorithm.
    auto cfg_fast = align_cfg::edit | align_cfg::aligned_ends{free_ends_first};

    // Computes semi global edit distance using slower standard pairwise algorithm.
    auto cfg_slow = align_cfg::edit | align_cfg::aligned_ends{free_ends_second};
//! [semi_global]
    (void) cfg_fast;
    (void) cfg_slow;
}

{
//! [max_error]
    using namespace seqan3;

    // Computes global edit distance allowing maximal 3 errors.
    auto cfg = align_cfg::edit | align_cfg::max_error{3u};
//! [max_error]
    (void) cfg;
}
}
