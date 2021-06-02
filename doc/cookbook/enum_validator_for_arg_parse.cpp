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

template <typename option_value_t>
class EnumValidator : public seqan3::value_list_validator<option_value_t>
{
public:
    //!\brief Type of values that are tested by validator
    using option_value_type = option_value_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */

    /*!\brief Constructing from a range.
     * \tparam range_type - the type of range; must model std::ranges::forward_range and
     *                      EnumValidator::option_value_type must be constructible from the rvalue reference type of the
     *                      given range
     * \param[in] rng - the range of valid values to test
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::constructible_from<option_value_type, std::ranges::range_rvalue_reference_t<range_type>>
    //!\endcond
    EnumValidator(range_type rng)
    {
        values.clear();
        std::ranges::move(std::move(rng), std::cpp20::back_inserter(values));
        std::ranges::sort(values);
        values.erase(std::unique(values.begin(), values.end()), values.end());
    }
    //!\}

    void operator()(option_value_type const & cmp) const {};

    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::convertible_to<std::ranges::range_value_t<range_type>, option_value_type>
    //!\endcond
    void operator()(range_type const & range) const {};

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        auto map = seqan3::enumeration_names<option_value_t>;
        std::vector<std::pair<std::string_view, option_value_t>> key_value_pairs(map.begin(), map.end());
        std::ranges::sort(key_value_pairs, [] (auto pair1, auto pair2)
            {
                if constexpr (std::totally_ordered<option_value_t>)
                {
                    if (pair1.second != pair2.second)
                        return pair1.second < pair2.second;
                }
                return pair1.first < pair2.first;
            });

        return seqan3::detail::to_string("Value must be one of (method name or number) ",
                                         key_value_pairs | std::views::keys,
                                         ".");
    }

private:

    //!\brief Minimum of the range to test.
    std::vector<option_value_type> values{};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    // Validatiors:
    EnumValidator<my_methods> method_validator{seqan3::enumeration_names<my_methods> | std::views::values};

    // Options:
    parser.add_option(args.methods, 'm', "method", "Choose the method(s) to be used. ",
                      seqan3::option_spec::standard, method_validator);
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
