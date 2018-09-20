//! [usage]
#include <seqan3/argument_parser/all.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("The-Age-App", argc, argv); // initialize

    int age{30}; // define default values directly in the variable

    myparser.add_option(age, 'a', "user-age", "Please specify your age.");

    try
    {
        myparser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
    {
        std::cerr << "The-Age-App - [PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        std::cout << std::endl << "Thanks for using The-Age-App!" << std::endl; // customize
        return 0;
    }

    std::cout << "integer given by user: " << age << std::endl;
    return 0;
}
//! [usage]
