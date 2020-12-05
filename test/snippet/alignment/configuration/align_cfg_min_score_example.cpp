#include <seqan3/alignment/configuration/align_config_min_score.hpp>
#include <seqan3/core/configuration/configuration.hpp>

int main()
{
    // Allow a minimal score of -5, i.e. at most 5 edit operations.
    seqan3::configuration config = seqan3::align_cfg::min_score{-5};
    auto min_score = std::get<seqan3::align_cfg::min_score>(config);
    min_score.score = -5;
}
