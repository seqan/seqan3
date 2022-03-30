#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

int main()
{
    // Give the seqan3 file type as a template argument to get all valid extensions for this file.
    seqan3::input_file_validator<seqan3::sequence_file_input<>> validator3{};
    seqan3::debug_stream << validator3.get_help_page_message() << '\n';

    return 0;
}
