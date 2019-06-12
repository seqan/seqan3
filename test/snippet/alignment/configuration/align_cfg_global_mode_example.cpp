#include <seqan3/alignment/configuration/align_config_mode.hpp>

int main()
{
    // Select the global mode.
    seqan3::align_cfg::mode cfg_global{seqan3::global_alignment};
}
