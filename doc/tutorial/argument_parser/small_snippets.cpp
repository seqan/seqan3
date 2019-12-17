#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/filesystem>

//![validator_include]
#include <seqan3/argument_parser/validators.hpp>
//![validator_include]

struct cmd_arguments
{
    std::filesystem::path file_path{};
    std::vector<uint8_t> seasons{};
    std::string aggregate_by{"mean"};
    bool header_is_set{};
};

cmd_arguments args{};

int main(int argc, char ** argv)
{

{
seqan3::argument_parser parser{"Example-Parser", argc, argv};
//![add_positional_option]
size_t variable{};
parser.add_positional_option(variable, "This is a description.");
//![add_positional_option]
}

{
seqan3::argument_parser parser{"Example-Parser", argc, argv};
//![add_option]
size_t variable{};
parser.add_option(variable, 'n', "my-number", "This is a description.");
//![add_option]
}

{
seqan3::argument_parser parser{"Example-Parser", argc, argv};
//![add_flag]
bool variable{};
parser.add_flag(variable, 'f', "my_flag", "This is a description.");
//![add_flag]
}

{
seqan3::argument_parser parser{"Example-Parser", argc, argv};
//![option_list]
std::vector<std::string> list_variable{};
parser.add_option(list_variable, 'n', "names", "Give me some names.");
//![option_list]
}

{
seqan3::argument_parser parser{"Example-Parser", argc, argv};
//![positional_option_list]
std::string variable{};
std::vector<std::string> list_variable{};
parser.add_positional_option(variable, "Give me a single variable.");
parser.add_positional_option(list_variable, "Give me one or more variables!.");
//![positional_option_list]
}

{
seqan3::argument_parser parser{"Example-Parser", argc, argv};
//![required_option]
std::string required_variable{};
parser.add_option(required_variable, 'n', "name", "I really need a name.", seqan3::option_spec::REQUIRED);
//![required_option]
}

{
seqan3::argument_parser parser{"Example-Parser", argc, argv};
//![input_file_validator]
parser.add_positional_option(args.file_path, "Please provide a tab separated data file.",
                             seqan3::input_file_validator{{"tsv"}});
//![input_file_validator]
}

}
