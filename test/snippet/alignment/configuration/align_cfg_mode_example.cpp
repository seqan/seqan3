#include <seqan3/alignment/configuration/align_config_mode.hpp>

int main()
{
//! [global]
    using namespace seqan3;

    // Select the global mode.
    align_cfg::mode cfg_global{align_cfg::global_alignment};
//! [global]

    (void) cfg_global;
}
