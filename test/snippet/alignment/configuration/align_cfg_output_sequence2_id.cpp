#include <seqan3/alignment/configuration/align_config_output.hpp>

int main()
{
    // Output only the id of the second sequence.
    seqan3::configuration cfg = seqan3::align_cfg::output_sequence2_id{};
}
