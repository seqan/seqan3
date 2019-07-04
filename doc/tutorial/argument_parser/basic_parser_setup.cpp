#include <seqan3/argument_parser/all.hpp>     // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>  // our custom output stream

using namespace seqan3;

int main(int argc, char ** argv)
{
    argument_parser myparser{"Game-of-Parsing", argc, argv};        // initialise myparser

    // ... add information, options, flags and positional options

    try
    {
         myparser.parse();                                          // trigger command line parsing
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "[Winter has come] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
}
