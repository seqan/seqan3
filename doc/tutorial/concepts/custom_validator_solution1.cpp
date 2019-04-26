#include <seqan3/argument_parser/all.hpp>

struct custom_validator
{
    using value_type = double; // used for all arithmetic types

    void operator() (value_type const &) const
    {
        // add implementation later
    }

    std::string get_help_page_message () const
    {
        // add real implementation later
        return "";
    }
};

static_assert(seqan3::Validator<custom_validator>); // does not cause compile error

int main() {}
