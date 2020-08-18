#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>

int main()
{
    // Computes semi global edit distance using fast-bit vector algorithm.
    auto cfg_fast = seqan3::align_cfg::method_global{
                        seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                        seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                        seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                        seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}} |
                    seqan3::align_cfg::edit_scheme |
                    seqan3::align_cfg::aligned_ends{seqan3::free_ends_first};

    // Computes semi global edit distance using slower standard pairwise algorithm.
    auto cfg_slow = seqan3::align_cfg::method_global{
                        seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                        seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                        seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                        seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}} |
                    seqan3::align_cfg::edit_scheme |
                    seqan3::align_cfg::aligned_ends{seqan3::free_ends_second};

    // Computes global edit distance allowing maximal 3 errors.
    auto cfg_errors = seqan3::align_cfg::method_global{} |
                      seqan3::align_cfg::edit_scheme |
                      seqan3::align_cfg::max_error{3u};
}
