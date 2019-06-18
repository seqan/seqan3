//! [argparse]
#include <string>                         // for std::string
#include <seqan3/argument_parser/all.hpp> // for argument_parser
#include <seqan3/core/debug_stream.hpp>   // for debug_stream

using namespace seqan3;

int main(int argc, char * argv[])
{
    // Create a buffer for the input.
    std::string input{};
    // Initialise the Argument Parser and add an option.
    seqan3::argument_parser parser("My first program", argc, argv);
    parser.add_positional_option(input, "Give me text.");

    // Parse the given arguments and catch possible errors.
    try
    {
        parser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext)
    {
        debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
        return 0;
    }

    debug_stream << "The text was: " << input << "\n";
}
//! [argparse]
