#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>

int main()
{
    // Configuration with linear gaps.
    seqan3::align_cfg::gap linear_cfg{seqan3::gap_scheme{seqan3::gap_score{-1}}};

    // Configuration with affine_gaps. Score for opening a gap during the alignment algorithm will be -11.
    seqan3::align_cfg::gap affine_cfg{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}}};
}
