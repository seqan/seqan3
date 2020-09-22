#include <seqan3/alignment/configuration/align_config_output.hpp>

int main()
{
    // Compute only the alignment.
    seqan3::configuration cfg = seqan3::align_cfg::output_alignment{};
}
