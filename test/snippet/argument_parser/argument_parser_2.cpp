#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"The-Age-App", argc, argv}; // initialize

    int age{30}; // define default values directly in the variable

    myparser.add_option(age, 'a', "user-age", "Please specify your age.");

    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "The-Age-App - [PARSER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    seqan3::debug_stream << "integer given by user: " << age << '\n';
    return 0;
}
