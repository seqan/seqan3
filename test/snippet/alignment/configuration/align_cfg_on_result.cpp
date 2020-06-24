#include <seqan3/alignment/configuration/align_config_on_result.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::align_cfg::on_result cfg{[] (auto && result) { seqan3::debug_stream << result << '\n'; }};
}
