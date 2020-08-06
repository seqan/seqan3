#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/configuration/on_result.hpp>

int main()
{
    seqan3::search_cfg::on_result cfg{[] (auto && result) { seqan3::debug_stream << result << '\n'; }};
    return 0;
}
