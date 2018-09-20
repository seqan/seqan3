//! [usage]
#include <seqan3/argument_parser/all.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Test", argc, argv); // initialize

    std::string myfile;
    seqan3::file_ext_validator my_validator{"fa","fasta"};

    myparser.add_option(myfile,'f',"file","Give me a filename.",
                        seqan3::option_spec::DEFAULT, my_validator);

    // an exception will be thrown if the user specifies a filename
    // that does not have one of the extensions ["fa","fasta"]
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

    std::cout << "filename given by user passed validation: " << myfile << "\n";
    return 0;
}
//! [usage]
