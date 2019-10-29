#include <seqan3/alignment/configuration/align_config_result.hpp>

int main()
{
    // Compute only the score.
    seqan3::align_cfg::result default_cfg{};
    // same as
    seqan3::align_cfg::result cfg_score{seqan3::with_score};

    // Compute the score and the back coordinate.
    seqan3::align_cfg::result cfg_end{seqan3::with_back_coordinate};

    // Compute the score, the back coordinate and the front coordinate.
    seqan3::align_cfg::result cfg_begin{seqan3::with_front_coordinate};

    // Compute the score, the back coordinate, the front coordinate and the alignment.
    seqan3::align_cfg::result cfg_alignment{seqan3::with_alignment};

    // You can also change the score type:

    // Compute only the score given a specific score_type.
    seqan3::align_cfg::result cfg_score_uint16{seqan3::with_score, seqan3::using_score_type<uint16_t>};

    // Compute the score given a specific score_type and the back coordinate.
    seqan3::align_cfg::result cfg_end_double{seqan3::with_back_coordinate, seqan3::using_score_type<double>};

    // ...
}
