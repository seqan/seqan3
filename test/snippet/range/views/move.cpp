#include <seqan3/core/platform.hpp>

#ifdef SEQAN3_DEPRECATED_310
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <string>

#include <seqan3/range/views/move.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

int main()
{
    std::vector<std::string> vec_in{"ABC", "DEF", "GEH"};

    std::vector<std::string> vec_out0{};
    vec_out0.resize(3);
    std::ranges::copy(vec_in,                       // copies strings from in to out
                      vec_out0.begin());

    std::vector<std::string> vec_out1{};
    vec_out1.resize(3);
    std::ranges::copy(vec_in | seqan3::views::move, // moves strings from in to out
                      vec_out1.begin());
}
#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310
