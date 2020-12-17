#include <seqan3/argument_parser/all.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv};
    std::string myvar{"Example"};
    myparser.add_option(myvar, 's', "special-op", "You know what you doin'?",
                        seqan3::option_spec::advanced);
}
