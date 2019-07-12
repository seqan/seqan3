#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

int main(int argc, char ** argv)
{
    argument_parser myparser{"Game-of-Parsing", argc, argv, false};
    // disable version checks permanently --------------------^
}
