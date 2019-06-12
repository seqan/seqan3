#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>

int main()
{
try
{
    // A symmetric band around the main diagonal.
    seqan3::align_cfg::band band_cfg{seqan3::static_band{seqan3::lower_bound{-4}, seqan3::upper_bound{4}}};

    // A band starting with the main diagonal shifted by 3 cells to the right.
    seqan3::align_cfg::band band_cfg_hi{seqan3::static_band{seqan3::lower_bound{3}, seqan3::upper_bound{7}}};

    // A band starting with the main diagonal shifted by 3 cells down.
    seqan3::align_cfg::band band_cfg_lo{seqan3::static_band{seqan3::lower_bound{-7}, seqan3::upper_bound{-3}}};

    // An invalid band configuration.
    seqan3::align_cfg::band band_cfg_invalid{seqan3::static_band{seqan3::lower_bound{7}, seqan3::upper_bound{3}}};
}
catch(...)
{
}
}
