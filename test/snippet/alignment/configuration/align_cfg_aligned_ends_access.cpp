#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    // Create an end_gaps object with one user defined static value and one user defined non-static value.
    seqan3::end_gaps eg{seqan3::front_end_first{std::true_type{}}, seqan3::front_end_second{true}};

    // Check if the front_end_first parameter contains static information.
    if constexpr (decltype(eg)::is_static<0>())
    {
        seqan3::debug_stream << "The leading gaps of the first sequence are static and the value is: " <<
                                std::boolalpha << decltype(eg)::get_static<0>() << '\n';
    }

    // Defaulted parameters will always be false and static.
    seqan3::debug_stream << "The leading gaps of the first sequence are static and the value is " <<
                            std::boolalpha << decltype(eg)::get_static<0>() << '\n';

    // Non-static parameters won't be captured as static. Trying to access it via get_static will trigger a static assert.
    if constexpr (!decltype(eg)::is_static<2>())
    {
        seqan3::debug_stream << "The leading gaps of the second sequence is not static! The value is: " <<
                                std::boolalpha << eg[2] << '\n';
    }

    seqan3::debug_stream << "The value can always be determined at runtime like for the trailing gaps of the second sequence: " <<
                            std::boolalpha << eg[3] << '\n';
}
