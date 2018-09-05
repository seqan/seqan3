//! [usage]
#include <seqan3/argument_parser/all.hpp>

int main(int argc, const char ** argv)
{
    seqan3::argument_parser myparser("Grade-Average", argc, argv); // initialize

    std::string name{"Max Muster"}; // define default values directly in the variable.
    bool bonus{false};
    std::vector<double> grades{};   // you can also specify a vector that is treated as a list option.

    myparser.add_option(name, 'n', "name", "Please specify your name.");
    myparser.add_flag(bonus, 'b', "bonus", "Please specify if you got the bonus.");
    myparser.add_positional_option(grades, "Please specify your grades.");

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

    if (bonus)
        grades.push_back(1.0); // extra good grade

    double avg{0};
    for (auto g : grades)
        avg += g;

    avg = avg / grades.size();

    std::cout << name << " has an average grade of " << avg << std::endl;

    return 0;
}
//! [usage]
