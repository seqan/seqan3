#include <system_error>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/ranges>

namespace seqan3::custom
{
// Specialise the seqan3::custom::argument_parsing data structure to enable parsing of std::errc.
template <>
struct argument_parsing<std::errc>
{
    // Specialise a mapping from an identifying string to the respective value of your type Foo.
    static inline std::unordered_map<std::string_view, std::errc> const enumeration_names
    {
            {"no_error", std::errc{}},
            {"timed_out", std::errc::timed_out},
            {"invalid_argument", std::errc::invalid_argument},
            {"io_error", std::errc::io_error}
    };
};

} // namespace seqan3::custom

int main(int argc, char const * argv[])
{
    std::errc value{};

    seqan3::argument_parser parser{"my_program", argc, argv};

    // Because of the argument_parsing struct and
    // the static member function enumeration_names
    // you can now add an option that takes a value of type std::errc:
    parser.add_option(value, 'e', "errc", "Give me a std::errc value.", seqan3::option_spec::DEFAULT,
                      seqan3::value_list_validator{(seqan3::enumeration_names<std::errc> | std::views::values)});

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR] " << ext.what() << "\n"; // customize your error message
        return -1;
    }

    return 0;
}
