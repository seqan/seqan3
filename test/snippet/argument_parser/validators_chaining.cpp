#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/filesystem>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv}; // initialize

    //![validator_call]
    std::string file_name;

    seqan3::regex_validator absolute_path_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"};
    seqan3::file_ext_validator my_file_ext_validator{{"sa", "so"}};

    myparser.add_option(file_name, 'f', "file","Give me a file name/path.",
                        seqan3::option_spec::DEFAULT, absolute_path_validator | my_file_ext_validator);
    //![validator_call]

    // an exception will be thrown if the user specifies a file name
    // that is not an absolute path or does not have one of the file extension [sa,so]
    try
    {
        myparser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    std::cout << "filename given by user passed validation: " << file_name << "\n";
    return 0;
}
