#include <seqan3/alignment/configuration/align_config_result.hpp>

int main()
{
//! [example]
    using namespace seqan3;

    // Compute only the score.
    align_cfg::result cfg_score{with_score};

    // Compute the score and the end coordinate.
    align_cfg::result cfg_end{with_back_coordinate};

    // Compute the score, the end coordinate and the begin coordinate.
    align_cfg::result cfg_begin{with_front_coordinate};

    // Compute the score, the end coordinate, the begin coordinate and the alignment.
    align_cfg::result cfg_alignment{with_alignment};

//! [example]
    (void) cfg_score;
    (void) cfg_end;
    (void) cfg_begin;
    (void) cfg_alignment;
}
