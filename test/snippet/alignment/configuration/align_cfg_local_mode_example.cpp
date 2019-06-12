#include <seqan3/alignment/configuration/align_config_mode.hpp>

int main()
{
    // Select the local mode.
    seqan3::align_cfg::mode cfg_local{seqan3::local_alignment};
}
