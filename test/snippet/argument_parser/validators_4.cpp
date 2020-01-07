#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv}; // initialize

    //![validator_call]
    std::string my_string;
    seqan3::regex_validator my_validator{"[a-zA-Z]+@[a-zA-Z]+\\.com"};

    myparser.add_option(my_string,'s',"str","Give me a string.",
                        seqan3::option_spec::DEFAULT, my_validator);
    //![validator_call]

    // an exception will be thrown if the user specifies a string
    // that is no email address ending on .com
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    seqan3::debug_stream << "email address given by user passed validation: " << my_string << "\n";
    return 0;
}
