#include <seqan3/std/algorithm>                         // for std::ranges::move
#include <unordered_map>                                // for std::unordered_map

#include <seqan3/argument_parser/argument_parser.hpp>   // for seqan3::argument_parser
#include <seqan3/argument_parser/auxiliary.hpp>         // for enumeration_names
#include <seqan3/argument_parser/exceptions.hpp>        // for seqan3::validation_error
#include <seqan3/core/debug_stream.hpp>                 // for seqan3::debug_stream

//!\brief An enum for the different methods.
enum my_methods
{
    method_a = 0,
    method_b = 1,
    method_c = 2,

    // Also add new methods to the default values in the argument parser

    //!\cond
    // ATTENTION: Must always be the last item; will be used to determine the number of ids.
    SIZE //!< Determines the size of the enum.
    //!\endcond
};

struct cmd_arguments
{
    std::vector<my_methods> methods{method_a, method_c};   // default: method a and c
};

std::unordered_map<std::string, my_methods> enumeration_names(my_methods)
{
return std::unordered_map<std::string, my_methods>{{"0", my_methods::method_a}, {"method_a", my_methods::method_a},
                                                   {"1", my_methods::method_b}, {"method_b", my_methods::method_b},
                                                   {"2", my_methods::method_c}, {"method_c", my_methods::method_c}};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    // Options:
    parser.add_option(args.methods, 'n', "nethod", "Choose the method(s) to be used. ",
                      seqan3::option_spec::standard);   // uses the default_validator
    parser.add_option(args.methods, 'm', "method", "Choose the method(s) to be used. ",
                      seqan3::option_spec::standard, seqan3::enum_validator{(seqan3::enumeration_names<my_methods>
                                                                             | std::views::values)});
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"myTool", argc, argv};             // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);

    // Parse the given arguments and catch possible errors.
    try
    {
        myparser.parse();                                               // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                   // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << '\n';       // customise your error message
        return -1;
    }

    // Check that method selection contains no duplicates.
    std::vector<my_methods> unique_methods{args.methods};
    std::ranges::sort(unique_methods);
    unique_methods.erase(std::unique(unique_methods.begin(), unique_methods.end()), unique_methods.end());
    if (args.methods.size() > unique_methods.size())
    {
        seqan3::debug_stream << "[Error] The same method was selected multiple times.\n";
        seqan3::debug_stream << "Methods to be used: " << args.methods << '\n';
        return -1;
    }

    return 0;
}
