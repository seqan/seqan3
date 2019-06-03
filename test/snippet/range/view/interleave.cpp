#include <range/v3/algorithm/equal.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/interleave.hpp>

int main()
{
    std::string u{"FOOBARBAXBAT"};
    std::string i{"in"};
    size_t s = 3;

    auto v = u | seqan3::view::interleave(s, i);

    seqan3::debug_stream << v << '\n'; // prints FOOinBARinBAXinBAT
}
