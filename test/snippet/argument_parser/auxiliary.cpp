#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

int main(int argc, const char ** argv)
{
//! [usage]
argument_parser myparser{"Test", argc, argv};
std::string myvar{"Example"};
myparser.add_option(myvar, 's', "special-op", "You know what you doin'?",
                    option_spec::ADVANCED);
//! [usage]
}
