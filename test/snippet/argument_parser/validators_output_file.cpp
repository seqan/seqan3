#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv}; // initialize

    //! [validator_call]
    std::filesystem::path myfile{};

    // Use the seqan3::output_file_open_options to indicate that you allow overwriting existing output files, ...
    myparser.add_option(myfile, 'f', "file", "Output file containing the processed sequences.",
                        seqan3::option_spec::standard,
                        seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"fa","fasta"}});

    // ... or that you will throw a seqan3::validation_error if the user specified output file already exists
    myparser.add_option(myfile, 'g', "file2", "Output file containing the processed sequences.",
                        seqan3::option_spec::standard,
                        seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"fa","fasta"}});
    //! [validator_call]

    // an exception will be thrown if the user specifies a filename
    // that does not have one of the extensions ["fa","fasta"],
    // if the file already exists, or if the file is not writable.
    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    seqan3::debug_stream << "filename given by user passed validation: " << myfile << "\n";
    return 0;
}
