#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>

int main()
{
//! [example]
    using namespace seqan3;

    // Configuration with linear gaps.
    align_cfg::gap linear_cfg{gap_scheme{gap_score{-1}}};

    // Configuration with affine_gaps. Score for opening a gap during the alignment algorithm will be -11.
    align_cfg::gap affine_cfg{gap_scheme{gap_score{-1}, gap_open_score{-10}}};
//! [example]

    (void) linear_cfg;
    (void) affine_cfg;
}
