#include <seqan3/argument_parser/all.hpp>
//![validator]
#include <cmath>

struct custom_validator
{
    using option_value_type = double; // used for all arithmetic types

    void operator() (option_value_type const & val) const
    {
        if ((std::round(val)                         != val) ||  // not an integer
            (std::pow(std::round(std::sqrt(val)), 2) != val))    // not a square
        {
            throw seqan3::validation_error{"The provided number is not an arithmetic square."};
        }
    }

    std::string get_help_page_message () const
    {
        return "Value must be the square of an integral number.";
    }
};
//![validator]

static_assert(seqan3::validator<custom_validator>);

//![main]
int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser("Test-Parser", argc, argv);

    int32_t variable{};
    int16_t variable2{};

    myparser.add_option(variable, 'i', "", "An int that is a square", seqan3::option_spec::DEFAULT,
                        custom_validator{}); // ← your validator is used!

    myparser.add_option(variable2, 'j', "", "An int that is a square and within [0,20].", seqan3::option_spec::DEFAULT,
                        custom_validator{} | seqan3::arithmetic_range_validator{0, 20}); // ← now it's chained

    try
    {
         myparser.parse(); // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << ext.what() << '\n';
        return -1;
    }

    seqan3::debug_stream << "Yeah!\n";

    return 0;
}
//![main]
