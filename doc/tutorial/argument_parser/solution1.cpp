#include <seqan3/argument_parser/all.hpp> // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>   // our custom output stream

using namespace seqan3;

void initialize_argument_parser(argument_parser & parser)
{
    parser.info.author = "Cersei";
    parser.info.short_description = "Aggregate average Game of Thrones viewers by season.";
    parser.info.version = "1.0.0";
}

int main(int argc, char ** argv)
{
    argument_parser myparser{"Game-of-Parsing", argc, argv};

    // code from assignment 1
}
