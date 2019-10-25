#include <seqan3/alignment/configuration/align_config_parallel.hpp>

int main()
{
    // Enables parallel computation with two threads.
    seqan3::align_cfg::parallel cfg{2};
}
