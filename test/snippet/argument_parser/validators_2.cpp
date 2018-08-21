//! [usage]
#include <seqan3/argument_parser/all.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Test", argc, argv); // initialize

    int myint;
    seqan3::value_list_validator my_validator{2,4,6,8,10};

    myparser.add_option(myint,'i',"integer","Give me a number.",
                        seqan3::option_spec::DEFAULT, my_validator);

    // an exception will be thrown if the user specifies an integer
    // that is not one of [2,4,6,8,10] (e.g. "./test_app -i 3")
    try
    {
        myparser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what(); // customize your error message
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        return 0;
    }

    std::cout << "integer given by user passed validation: " << myint << "\n";
    return 0;
}
//! [usage]
