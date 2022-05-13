#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>     // our custom output stream

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"Game-of-Parsing", argc, argv};        // initialise myparser

    // ... add information, options, flags and positional options

    try
    {
         myparser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "[Winter has come] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
}
