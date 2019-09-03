#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

int main()
{
    // Here we use the global and banded alignment configurations to show how they can be combined.
    seqan3::configuration my_cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
                                   seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-4},
                                                                               seqan3::upper_bound{4}}};
}
