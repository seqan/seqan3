#include <seqan3/alignment/configuration/align_config_mode.hpp>

int main()
{
//! [global]
    using namespace seqan3;

    // Select the global mode.
    align_cfg::mode cfg_global{global_alignment};
//! [global]

//! [local]
    // Select the local mode.
    align_cfg::mode cfg_local{local_alignment};
//! [local]

    (void) cfg_global;
    (void) cfg_local;
}
