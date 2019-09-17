#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <range/v3/view/zip_with.hpp>


int main()
{
    std::vector<std::string> id{};

    std::vector<std::string> src{"hello", "world"};

    auto v = ranges::view::zip_with([&](auto cur){ id.push_back(std::move(cur)); }, id);
    std::ranges::copy(src, std::ranges::begin(v));

    seqan3::debug_stream << id << '\n';
}
