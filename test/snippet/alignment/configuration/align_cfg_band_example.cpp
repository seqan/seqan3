#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>

int main()
{
try
{
//! [example]
    using namespace seqan3;

    // A symmetric band around the main diagonal.
    align_cfg::band band_cfg{static_band{lower_bound{-4}, upper_bound{4}}};

    // A band starting with the main diagonal shifted by 3 cells to the right.
    align_cfg::band band_cfg_hi{static_band{lower_bound{3}, upper_bound{7}}};

    // A band starting with the main diagonal shifted by 3 cells down.
    align_cfg::band band_cfg_lo{static_band{lower_bound{-7}, upper_bound{-3}}};

    // An invalid band configuration.
    align_cfg::band band_cfg_invalid{static_band{lower_bound{7}, upper_bound{3}}};
//! [example]

    (void) band_cfg;
    (void) band_cfg_hi;
    (void) band_cfg_lo;
    (void) band_cfg_invalid;
}
catch(...)
{
}
}
