//! [usage]
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Test", argc, argv); // initialize

    std::string my_string;
    seqan3::regex_validator my_validator{"[a-zA-Z]+@[a-zA-Z]+\\.com"};

    myparser.add_option(my_string,'s',"str","Give me a string.",
                        seqan3::option_spec::DEFAULT, my_validator);

    // an exception will be thrown if the user specifies a string
    // that is no email address ending on .com
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // catch all other exceptions caused by user errors
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        return 0;
    }

    seqan3::debug_stream << "email address given by user passed validation: " << my_string << "\n";
    return 0;
}
//! [usage]
