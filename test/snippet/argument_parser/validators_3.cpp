//! [usage]
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

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
    catch (seqan3::argument_parser_error const & ext) // catch all other exceptions caused by user errors
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        return 0;
    }

    seqan3::debug_stream << "filename given by user passed validation: " << myfile << "\n";
    return 0;
}
//! [usage]
