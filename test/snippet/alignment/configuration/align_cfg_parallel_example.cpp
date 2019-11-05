#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <thread>

int main()
{
    // Enables parallel computation with two threads.
    seqan3::align_cfg::parallel cfg_2{2};

    // Enables parallel computation with the number of concurrent threads supported by the current architecture.
    seqan3::align_cfg::parallel cfg_n{std::thread::hardware_concurrency()};
    ;
}
