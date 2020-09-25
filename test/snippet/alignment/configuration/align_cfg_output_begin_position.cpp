#include <seqan3/alignment/configuration/align_config_output.hpp>

int main()
{
    // Compute only the begin position of the aligned sequences.
    seqan3::configuration cfg = seqan3::align_cfg::output_begin_position{};
}
