//! [usage]
#include <seqan3/argument_parser/all.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Penguin_Parade", argc, argv); // initialize

    myparser.info.version = "2.0.0";
    myparser.info.date = "12.01.2017";
    myparser.info.short_description = "Organize your penguin parade";
    myparser.info.description.push_back("First Paragraph.");
    myparser.info.description.push_back("Second Paragraph.");
    myparser.info.examples.push_back("./penguin_parade Skipper Kowalski Rico Private -d 10 -m 02 -y 2017");

    int d{01};   // day
    int m{01};   // month
    int y{2050}; // year

    myparser.add_option(d, 'd', "day", "Please specify your preferred day.");
    myparser.add_option(m, 'm', "month", "Please specify your preferred month.");
    myparser.add_option(y, 'y', "year", "Please specify your preferred year.");

    std::vector<std::string> penguin_names;

    myparser.add_positional_option(penguin_names, "Specify the names of the penguins.");

    try
    {
        myparser.parse();
    }
    catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
    {
        std::cerr << ext.what();
        return -1;
    }
    catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
    {
        return 0;
    }

    // organize ...

    return 0;
}
//! [usage]
