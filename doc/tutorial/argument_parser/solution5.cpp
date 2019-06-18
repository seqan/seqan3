#include <fstream>
#include <numeric>

#include <range/v3/view/split.hpp>

#include <seqan3/argument_parser/all.hpp> // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>   // our custom output stream
#include <seqan3/std/charconv>            // includes std::from_chars
#include <seqan3/std/filesystem>          // use std::filesystem::path

using namespace seqan3;

// This is the program!
// Take a look at it if you are interested in an example of parsing a data file.
// -----------------------------------------------------------------------------
template <typename number_type, typename range_type>
number_type to_number(range_type && range)
{
    std::string str;
    number_type num;
    std::ranges::copy(range, std::back_inserter(str));
    auto res = std::from_chars(&str[0], &str[0] + str.size(), num);

    if (res.ec != std::errc{})
    {
        debug_stream << "Could not cast '" << range << "' to a valid number\n";
        throw std::invalid_argument{"CAST ERROR"};
    }
    return num;
}

void run_program(std::filesystem::path & path, std::vector<uint8_t> sn, std::string & aggr_by, bool hd_is_set)
{
    std::ifstream file{path.string()};

    if (file.is_open())
    {
        std::vector<double> v;
        std::string line;

        if (hd_is_set)
            std::getline(file, line); // ignore first line

        while (std::getline(file, line))
        {
            auto splitted_line = line | std::view::split('\t');
            auto it = splitted_line.begin(); // move to 1rst column

            if (std::find(sn.begin(), sn.end(), to_number<uint8_t>(*it)) != sn.end())
                v.push_back(to_number<double>(*std::next(it, 4)));
        }

        if (aggr_by == "median")
            debug_stream << ([&v] () { std::sort(v.begin(), v.end()); return v[v.size()/2]; })() << std::endl;
        else if (aggr_by == "mean")
            debug_stream << ([&v] () { double sum{}; for (auto i : v) sum += i; return sum / v.size(); })() << std::endl;
        else
            debug_stream << "I do not know the aggregation method " << aggr_by << std::endl;
    }
    else
    {
        debug_stream << "Error: Cannot open file for reading.\n";
    }
}
// -----------------------------------------------------------------------------

struct cmd_arguments
{
    std::filesystem::path file_path{};
    std::vector<uint8_t> seasons{};
    std::string aggregate_by{"mean"};
    bool header_is_set{};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Cercei";
    parser.info.short_description = "Aggregate average Game of Thrones viewers by season.";
    parser.info.version = "1.0.0";

    parser.add_positional_option(args.file_path, "Please provide a tab separated data file.");

//![solution]
    parser.add_option(args.seasons, 's', "season", "Choose the seasons to aggregate.", option_spec::REQUIRED);
//![solution]

    parser.add_option(args.aggregate_by, 'a', "aggregate-by", "Choose your method of aggregation: mean or median.");

    parser.add_flag(args.header_is_set, 'H', "header-is-set", "Let us know whether your data file contains a "
                                                              "header to ensure correct parsing.");
}

int main(int argc, char ** argv)
{
    argument_parser myparser{"Game of Parsing", argc, argv};        // initialise myparser
    cmd_arguments args{};

    initialize_argument_parser(myparser, args);

    try
    {
         myparser.parse();                                          // trigger command line parsing
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "[Winter has come] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    // parsing was successful !
    // we can start running our program
    run_program(args.file_path, args.seasons, args.aggregate_by, args.header_is_set);

    return 0;
}
