#include <seqan3/alignment/configuration/align_config_band.hpp>

int main()
{
try
{
    // A symmetric band around the main diagonal.
    seqan3::align_cfg::band_fixed_size band_cfg{seqan3::align_cfg::lower_diagonal{-4},
                                                seqan3::align_cfg::upper_diagonal{4}};

    // A band starting with the main diagonal shifted by 3 cells to the right.
    seqan3::align_cfg::band_fixed_size band_cfg_hi{seqan3::align_cfg::lower_diagonal{3},
                                                   seqan3::align_cfg::upper_diagonal{7}};

    // A band starting with the main diagonal shifted by 3 cells down.
    seqan3::align_cfg::band_fixed_size band_cfg_lo{seqan3::align_cfg::lower_diagonal{-7},
                                                   seqan3::align_cfg::upper_diagonal{-3}};

    // An invalid band configuration.
    seqan3::align_cfg::band_fixed_size band_cfg_invalid{seqan3::align_cfg::lower_diagonal{7},
                                                        seqan3::align_cfg::upper_diagonal{3}};
}
catch(...)
{
}
}
