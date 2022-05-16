#include <seqan3/argument_parser/all.hpp>

struct custom_validator
{
    using option_value_type = double; // used for all arithmetic types

    void operator() (option_value_type const &) const
    {
        // add implementation later
    }

    std::string get_help_page_message () const
    {
        // add real implementation later
        return "";
    }
};

static_assert(seqan3::validator<custom_validator>); // does not cause compile error

int main() {}
