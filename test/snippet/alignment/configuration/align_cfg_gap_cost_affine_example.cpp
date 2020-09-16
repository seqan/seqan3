#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>

int main()
{
    // Configuration with linear gap costs.
    seqan3::align_cfg::gap_cost_affine linear_cfg{seqan3::align_cfg::open_score{0},
                                                  seqan3::align_cfg::extension_score{-1}};

    // Configuration with affine gap costs. Score for opening a gap during the alignment algorithm will be -11.
    seqan3::align_cfg::gap_cost_affine affine_cfg{seqan3::align_cfg::open_score{-1},
                                                  seqan3::align_cfg::extension_score{-10}};

    // Accessing the members of the gap scheme
    int open = affine_cfg.open_score; // open == -1
    int extension = affine_cfg.extension_score; // extension == -10
}
