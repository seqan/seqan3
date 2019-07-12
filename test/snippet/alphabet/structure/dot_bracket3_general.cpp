#include <vector>

#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_db3;

    // create vector
    std::vector<seqan3::dot_bracket3> vec{'.'_db3, ')'_db3, ')'_db3};

    // modify and print
    vec[1] = '('_db3;

    for (seqan3::dot_bracket3 chr : vec)
        seqan3::debug_stream << seqan3::to_char(chr);  // .()

    seqan3::debug_stream << "\n";
}
