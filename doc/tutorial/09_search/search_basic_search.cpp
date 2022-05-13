#include <seqan3/core/debug_stream.hpp>       // pretty printing
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

using namespace std::string_literals; // for using the ""s string literal

int main()
{
    std::string text{"Garfield the fat cat without a hat."};
    seqan3::fm_index index{text};
    seqan3::debug_stream << search("cat"s, index) << '\n'; // [<query_id:0, reference_id:0, reference_pos:17>]
}
