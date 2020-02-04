#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"Grade-Average", argc, argv}; // initialize

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
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    if (bonus)
        grades.push_back(1.0); // extra good grade

    double avg{0};
    for (auto g : grades)
        avg += g;

    avg = avg / grades.size();

    seqan3::debug_stream << name << " has an average grade of " << avg << '\n';

    return 0;
}
