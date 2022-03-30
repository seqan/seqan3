#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    // Give the seqan3 file type as a template argument to get all valid extensions for this file.
    seqan3::output_file_validator<seqan3::sequence_file_output<>> validator3{
        seqan3::output_file_open_options::create_new};
    seqan3::debug_stream << validator3.get_help_page_message() << '\n';

    return 0;
}
