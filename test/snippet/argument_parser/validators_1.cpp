//! [usage]
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Test", argc, argv); // initialize

    int myint;
    seqan3::integral_range_validator<int> my_validator{2, 10};

    myparser.add_option(myint,'i',"integer","Give me a number.",
                        seqan3::option_spec::DEFAULT, my_validator);

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
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        return 0;
    }

    seqan3::debug_stream << "integer given by user passed validation: " << myint << "\n";
    return 0;
}
//! [usage]
