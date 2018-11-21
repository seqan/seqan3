//! [usage]
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser("The-Age-App", argc, argv); // initialize

    int age{30}; // define default values directly in the variable

    myparser.add_option(age, 'a', "user-age", "Please specify your age.");

    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // catch all other exceptions caused by user errors
    {
        std::cerr << "The-Age-App - [PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        seqan3::debug_stream << std::endl << "Thanks for using The-Age-App!" << std::endl; // customize
        return 0;
    }

    seqan3::debug_stream << "integer given by user: " << age << std::endl;
    return 0;
}
//! [usage]
