#include <seqan3/core/char_operations/predicate.hpp>// for is_char
#include <seqan3/core/debug_stream.hpp>             // for debug_stream
#include <seqan3/range/view/single_pass_input.hpp>  // for view::single_pass_input
#include <seqan3/range/view/take_until.hpp>         // for view::take_until*
#include <seqan3/std/ranges>                        // for std::view::reverse

using namespace seqan3;

int main()
{
    // regular usage
    std::string vec{"foo\nbar"};
    auto v = vec | view::take_until(is_char<'\n'>); // or use a lambda
    debug_stream << v << '\n'; // "foo"

    auto v2 = vec | std::view::reverse | view::take_until(is_char<'\n'>);
    debug_stream << v2 << '\n'; // "rab"

    // consuming behaviour
    std::string vec2{"foo      bar"}; // ← multiple spaces
    auto vin = vec2 | view::single_pass_input;
    auto v3 = vin | view::take_until_and_consume(is_blank);
    debug_stream << v3 << '\n'; // "foo"
    debug_stream << *std::ranges::begin(vin) << '\n'; // "b", the spaces where skipped
}
