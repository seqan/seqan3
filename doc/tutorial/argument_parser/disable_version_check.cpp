#include <seqan3/argument_parser/all.hpp>

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"Game-of-Parsing", argc, argv, false};
    // disable version checks permanently ----------------------------^
}
