//! [usage]
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv}; // initialize

    //![validator_call]
    int myint;
    seqan3::arithmetic_range_validator my_validator{2, 10};

    myparser.add_option(myint,'i',"integer","Give me a number.",
                        seqan3::option_spec::DEFAULT, my_validator);
    //![validator_call]

    // an exception will be thrown if the user specifies an integer
    // that is not in range [2,10] (e.g. "./test_app -i 15")
    try
    {
        myparser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    seqan3::debug_stream << "integer given by user passed validation: " << myint << "\n";
    return 0;
}
//! [usage]
