//! [usage]
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Test", argc, argv); // initialize

    seqan3::filesystem::path myfile;
    seqan3::file_existance_validator my_validator{};

    myparser.add_option(myfile,'f',"file","Give me a filename.",
                        seqan3::option_spec::DEFAULT, my_validator);

    // an exception will be thrown if the user specifies a filename that does not exist
    try
    {
        myparser.parse();
    }
    catch (seqan3::validation_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        return 0;
    }

    seqan3::debug_stream << "filename given by user passed validation: " << myfile.string() << "\n";

    return 0;
}
//! [usage]
