#include <string>

#include <seqan3/core/char_operations/predicate.hpp>// for is_char
#include <seqan3/core/debug_stream.hpp>             // for debug_stream
#include <seqan3/range/view/single_pass_input.hpp>  // for view::single_pass_input
#include <seqan3/range/view/take_until.hpp>         // for view::take_until*
#include <seqan3/std/ranges>                        // for std::view::reverse

int main()
{
    // regular usage
    std::string vec{"foo\nbar"};
    auto v = vec | seqan3::view::take_until(seqan3::is_char<'\n'>); // or use a lambda
    seqan3::debug_stream << v << '\n'; // "foo"

    auto v2 = vec | std::view::reverse | seqan3::view::take_until(seqan3::is_char<'\n'>);
    seqan3::debug_stream << v2 << '\n'; // "rab"

    // consuming behaviour
    std::string vec2{"foo      bar"}; // ← multiple spaces
    auto vin = vec2 | seqan3::view::single_pass_input;
    auto v3 = vin | seqan3::view::take_until_and_consume(seqan3::is_blank);
    seqan3::debug_stream << v3 << '\n'; // "foo"
    seqan3::debug_stream << *std::ranges::begin(vin) << '\n'; // "b", the spaces where skipped
}
