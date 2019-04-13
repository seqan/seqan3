#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/std/filesystem>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Test", argc, argv); // initialize

    std::filesystem::path myfile;

    myparser.add_option(myfile,'d',"dir","The output directory for storing the files.",
                        seqan3::option_spec::DEFAULT, seqan3::output_directory_validator());

    // an exception will be thrown if the user specifies a filename
    // that does not have one of the extensions ["fa","fasta"]
    try
    {
        myparser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    seqan3::debug_stream << "filename given by user passed validation: " << myfile << "\n";
    return 0;
}
