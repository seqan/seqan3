#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/tuple/pod_tuple.hpp>

int main()
{
    seqan3::pod_tuple<int, float> t{3, 4.7};
    static_assert(std::is_standard_layout_v<seqan3::pod_tuple<int, float>>);
    static_assert(std::is_trivial_v<seqan3::pod_tuple<int, float>>);

    // template parameters are automatically deduced:
    seqan3::pod_tuple t2{17, 3.7f, 19l};

    seqan3::debug_stream << std::get<0>(t2) << '\n'; // 17

    auto [ i, f, l ] = t2; // creates an int i with value 17, float f...
}
