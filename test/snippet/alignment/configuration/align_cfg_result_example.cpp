#include <seqan3/alignment/configuration/align_config_result.hpp>

int main()
{
//! [example]
    using namespace seqan3;

    // Compute only the score.
    align_cfg::result cfg_score{with_score};

    // Compute the score and the back coordinate.
    align_cfg::result cfg_end{with_back_coordinate};

    // Compute the score, the back coordinate and the front coordinate.
    align_cfg::result cfg_begin{with_front_coordinate};

    // Compute the score, the back coordinate, the front coordinate and the alignment.
    align_cfg::result cfg_alignment{with_alignment};

//! [example]
    (void) cfg_score;
    (void) cfg_end;
    (void) cfg_begin;
    (void) cfg_alignment;
}
