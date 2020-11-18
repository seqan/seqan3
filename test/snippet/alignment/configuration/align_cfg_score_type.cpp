#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/core/configuration/configuration.hpp>

int main()
{
    // Compute only the score.
    seqan3::configuration cfg1 = seqan3::align_cfg::score_type<int16_t>{}; // Now the alignment computes 16 bit integers.
    seqan3::configuration cfg2 = seqan3::align_cfg::score_type<float>{};  // Now the alignment computes float scores.
}
