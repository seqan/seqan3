#include <seqan3/alignment/configuration/align_config_max_error.hpp>

int main()
{
    // Allow maximal 5 errors.
    seqan3::align_cfg::max_error cfg{5u};
}
