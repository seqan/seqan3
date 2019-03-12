#include <seqan3/alignment/configuration/align_config_max_error.hpp>

int main()
{
//! [example]
    using namespace seqan3;

    // Allow maximal 5 errors.
    align_cfg::max_error cfg{5u};
//! [example]

    (void) cfg;
}
