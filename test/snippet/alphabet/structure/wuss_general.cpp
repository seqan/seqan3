#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_wuss51;

    // create vector
    std::vector<seqan3::wuss51> vec{'.'_wuss51, '>'_wuss51, '>'_wuss51};
    // modify and print
    vec[1] = '<'_wuss51;
    for (seqan3::wuss51 chr : vec)
        seqan3::debug_stream << seqan3::to_char(chr);  // .<>
    seqan3::debug_stream << "\n";
}
