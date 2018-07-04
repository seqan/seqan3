// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Contains some standard validators for (positional) options.
 */

#pragma once

#include <sstream>
#include <regex>

#include <range/v3/view/all.hpp> // ranges::view::all()

#include <seqan3/argument_parser/auxiliary.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/concept/core_language.hpp>
#include <seqan3/std/concept/callable.hpp>

namespace seqan3
{

/*!\interface seqan3::validator_concept <>
 * \brief The concept for option validators passed to add_option/positional_option.
 * \ingroup argument_parser
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename validator_type>
concept bool validator_concept = invocable_concept<validator_type, typename validator_type::value_type> &&
                                 requires(validator_type validator,
                                          typename validator_type::value_type value)
{
    typename validator_type::value_type;

    { validator.get_help_page_message() } -> std::string;
};
//!\endcond

//!\cond
template <typename option_value_type>
class integral_range_validator;
//!\endcond

/*!\brief A validator that checks whether a number is inside a given range.
 * \ingroup argument_parser
 *
 * \tparam option_value_type Must be a (container of) integral type(s).
 *
 * \details
 *
 * On construction, the validator must receive a maximum and a minimum number.
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever a given value does not lie inside the given min/max range.
 *
 *```cpp
 *     int main(int argc, char ** argv)
 *     {
 *          seqan3::argument_parser myparser(argc, argv); // initialize
 *
 *          int myint;
 *          seqan3::integral_range_validator my_validator<int>{2, 10};
 *
 *          myparser.add_option(myint,'i',"integer","Give me a number.",
 *                              seqan3::option_spec::DEFAULT, my_validator);
 *
 *          // an exception will be thrown if the user specifies an integer
 *          // that is not in range [2,10] (e.g. "./test_app -i 15")
 *          try
 *          {
 *               myparser.parse();
 *          }
 *          catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
 *          {
 *               std::cerr << "[PARSER ERROR] " << ext.what(); // customize your error message
 *               return -1;
 *          }
 *          catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
 *          {
 *               return 0;
 *          }
 *
 *          std::cout << "integer given by user passed validation: " << myint << endl;
 *          return 0;
 *     }
 *```
 */
template <integral_concept option_value_type>
class integral_range_validator<option_value_type>
{
public:
    //!\brief The type of value that this validator invoked upon.
    using value_type = option_value_type;

    /*!\brief The constructor.
     * \param[in] min_ Minimum set for the range to test.
     * \param[in] max_ Maximum set for the range to test.
     */
    integral_range_validator(value_type const min_, value_type const max_) :
        min{min_}, max{max_}
    {}

    /*!\brief Tests whether cmp lies inside [min,max].
     * \param cmp The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(value_type const & cmp) const
    {
        if (!((cmp <= max) && (cmp >= min)))
            throw parser_invalid_argument(detail::to_string("Validation Failed - Value ", cmp,
                                                            " is not in range [", min, ",", max, "]."));
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must be in range [", min, ",", max, "].");
    }

private:
    //!\brief Minimum of the range to test.
    value_type min{};

    //!\brief Maximum of the range to test.
    value_type max{};
};

template <container_concept option_value_type>
//!\cond
    requires integral_concept<typename option_value_type::value_type>
//!\endcond
class integral_range_validator<option_value_type>
{
public:
    //!\brief Type of values that are tested by validator (container)
    using value_type = option_value_type;
    //!\brief Underlying type of the container
    using inner_value_type = typename value_type::value_type;

    /*!\brief The constructor.
     * \param[in] min_ Minimum set for the range to test.
     * \param[in] max_ Maximum set for the range to test.
     */
    integral_range_validator(inner_value_type const min_,
                             inner_value_type const max_) :
        min{min_}, max{max_}
    {}

    /*!\brief Tests whether cmp lies inside [min,max].
     * \param cmp The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(value_type const & cmp) const
    {
        std::for_each(cmp.begin(), cmp.end(), [&] (auto cmp_v)
            {
                if (!((cmp_v <= max) && (cmp_v >= min)))
                    throw parser_invalid_argument(detail::to_string("Validation Failed - Value ", cmp_v,
                                                                    " is not in range [", min, ",", max, "]."));
            });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must be in range [", min, ",", max, "].");
    }

private:
    //!\brief Minimum of the range to test.
    inner_value_type min{};

    //!\brief Maximum of the range to test.
    inner_value_type max{};
};

/*!\brief A validator that checks whether a value is inside a list of valid values.
 * \ingroup argument_parser
 *
 * \details
 *
 * On construction, the validator must receive a list (vector) of valid values.
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever a given value is not in the given list.
 *
 *```cpp
 *     int main(int argc, char ** argv)
 *     {
 *          seqan3::argument_parser myparser(argc, argv); // initialize
 *
 *          int myint;
 *          seqan3::value_list_validator my_validator{2,4,6,8,10};
 *
 *          myparser.add_option(myint,'i',"integer","Give me a number.",
 *                              seqan3::option_spec::DEFAULT, my_validator);
 *
 *          // an exception will be thrown if the user specifies an integer
 *          // that is not one of [2,4,6,8,10] (e.g. "./test_app -i 3")
 *          try
 *          {
 *               myparser.parse();
 *          }
 *          catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
 *          {
 *               std::cerr << "[PARSER ERROR] " << ext.what(); // customize your error message
 *               return -1;
 *          }
 *          catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
 *          {
 *               return 0;
 *          }
 *
 *          std::cout << "integer given by user passed validation: " << myint << endl;
 *          return 0;
 *     }
 *```
 */
template <typename option_value_type>
class value_list_validator
{
public:
    //!\brief Type of values that are tested by validator
    using value_type = option_value_type;

    /*!\brief Constructing from a vector.
     * \param[in] v The vector of valid values to test.
     */
    value_list_validator(std::vector<option_value_type> const & v) :
        values{v}
    {}

    /*!\brief Constructing from an initializer_list.
     * \param[in] v The initializer_list of valid values to test.
     */
    value_list_validator(std::initializer_list<option_value_type> const & v) :
        values{v}
    {}

    /*!\brief Tests whether cmp lies inside values.
     * \param cmp The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(option_value_type const & cmp) const
    {
        if (!(std::find(values.begin(), values.end(), cmp) != values.end()))
            throw parser_invalid_argument(detail::to_string("Value ", cmp, " is not one of ", ranges::view::all(values), "."));
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must be one of ", ranges::view::all(values), ".");
    }

private:

    //!\brief Minimum of the range to test.
    std::vector<option_value_type> values;
};

/*!\brief A validator that checks if each value in a container appears in a list of valid values.
 * \ingroup argument_parser
 * \extends seqan3::value_list_validator
 *
 * \tparam option_value_type The container type. Must satisfy the seqan3::container_concept.
 */
template <container_concept option_value_type>
//!\cond
    requires !std::is_same_v<option_value_type, std::string>
//!\endcond
class value_list_validator<option_value_type>
{
public:
    //!\brief Type of values that are tested by validator (container)
    using value_type = option_value_type;
    //!\brief Underlying type of the container
    using inner_value_type = typename value_type::value_type;

    /*!\brief Constructing from a vector.
     * \param[in] v The vector of valid values to test.
     */
    value_list_validator(std::vector<inner_value_type> const & v) :
        values{v}
    {}

    /*!\brief Constructing from an initializer_list.
     * \param[in] v The initializer_list of valid values to test.
     */
    value_list_validator(std::initializer_list<inner_value_type> const & v) :
        values{v}
    {}

    /*!\brief Tests whether cmp lies inside values.
     * \param cmp The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(value_type const & cmp) const
    {
        std::for_each(cmp.begin(), cmp.end(), [&] (auto cmp_v)
            {
                if (!(std::find(values.begin(), values.end(), cmp_v) != values.end()))
                    throw parser_invalid_argument(detail::to_string("Value ", cmp_v, " is not one of ",
                                                                ranges::view::all(values), "."));
            });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must be one of ", ranges::view::all(values), ".");
    }

private:

    //!\brief Minimum of the range to test.
    std::vector<inner_value_type> values;
};

/*!\brief A validator that checks if a filenames has one of the valid extensions.
 * \ingroup argument_parser
 *
 * \details
 *
 * On construction, the validator must receive a list (vector) of valid file extensions.
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever a given filename (string) is not in the given extension list.
 *
 *```cpp
 * int main(int argc, char ** argv)
 * {
 *      seqan3::argument_parser myparser(argc, argv); // initialize
 *
 *      std::string myfile;
 *      seqan3::file_ext_validator my_validator{"fa","fasta"};
 *
 *      myparser.add_option(myfile,'f',"file","Give me a filename.",
 *                          seqan3::option_spec::DEFAULT, my_validator);
 *
 *      // an exception will be thrown if the user specifies a filename
 *      // that does not have one of the extensions ["fa","fasta"]
 *      try
 *      {
 *           myparser.parse();
 *      }
 *      catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
 *      {
 *           std::cerr << "[PARSER ERROR] " << ext.what(); // customize your error message
 *           return -1;
 *      }
 *      catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
 *      {
 *           return 0;
 *      }
 *
 *      std::cout << "filename given by user passed validation: " << myfile << endl;
 *      return 0;
 * }
 *```
 */
class file_ext_validator // TODO AFTER INITIAL MERGE: Check if file exists, also allow filepath
{
public:
    //!\brief Type of values that are tested by validator
    using value_type = std::string;

    /*!\brief Constructing from a vector.
     * \param[in] v The vector of valid file extensions to test (e.g. {"fa", "fasta"}).
     */
    file_ext_validator(std::vector<std::string> const & v) :
        extensions{v}
    {}

    /*!\brief Constructing from an initializer_list.
     * \param[in] v The initializer_list of valid file extensions to test (e.g. {"fa", "fasta"}).
     */
    file_ext_validator(std::initializer_list<std::string> const & v) :
        extensions{v}
    {}

    /*!\brief Tests whether cmp lies inside extensions.
     * \param cmp The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(std::string const & cmp) const
    {
        std::string ext{cmp.substr(cmp.find_last_of('.') + 1)};

        if (!(std::find(extensions.begin(), extensions.end(), ext) != extensions.end()))
            throw parser_invalid_argument(detail::to_string("Extension ", ext, " is not one of ",
                                                            ranges::view::all(extensions), "."));
    }

    //!\brief Tests whether every value of v lies inside extensions.
    void operator()(std::vector<std::string> const & v) const
    {
        std::for_each(v.begin(), v.end(), [&] (auto cmp) { (*this)(cmp); });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("File name extension must be one of ", ranges::view::all(extensions), ".");
    }

private:
    //!\brief Stores valid file extensions.
    std::vector<std::string> extensions;
};

//!\cond
template <typename option_value_type>
class regex_validator;
//!\endcond

/*!\brief A validator that checks if a matches a regular expression pattern.
 * \ingroup argument_parser
 *
 * \details
 *
 * On construction, the validator must receive a pattern for a regular expression.
 * The pattern variable will be used for constructing an std::regex and the
 * validator will call std::regex_match on the command line argument.
 * Note: A regex_match will only return true if the strings matches the pattern
 * completely (in contrast to regex_search which also matches substrings).
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever a given filename (string) is not in the given extension list.
 *
 *```cpp
 * int main(int argc, char ** argv)
 * {
 *      seqan3::argument_parser myparser(argc, argv); // initialize
 *
 *      std::string my_string;
 *      seqan3::regex_validator my_validator{"[a-zA-Z]+@[a-zA-Z]+\\.com"};
 *
 *      myparser.add_option(my_string,'s',"str","Give me a string.",
 *                          seqan3::option_spec::DEFAULT, my_validator);
 *
 *      // an exception will be thrown if the user specifies a string
 *      // that is no email address ending on .com
 *      try
 *      {
 *           myparser.parse();
 *      }
 *      catch (seqan3::parser_invalid_argument const & ext) // the user did something wrong
 *      {
 *           std::cerr << "[PARSER ERROR] " << ext.what(); // customize your error message
 *           return -1;
 *      }
 *      catch (seqan3::parser_interruption const &) // expected behaviour on special requests (e.g. `--help`)
 *      {
 *           return 0;
 *      }
 *
 *      std::cout << "email address given by user passed validation: " << my_string << endl;
 *      return 0;
 * }
 *```
 */
template <>
class regex_validator<std::string>
{
public:
    //!\brief Type of values that are tested by validator.
    using value_type = std::string;

    /*!\brief Constructing from a vector.
     * \param[in] pattern_ The pattern to match.
     */
    regex_validator(std::string const & pattern_) :
        pattern{pattern_}
    {}

    /*!\brief Tests whether cmp lies inside values.
     * \param[in] cmp The value to validate.
     * \throws parser_invalid_argument
     */
    void operator()(value_type const & cmp) const
    {
        std::regex rgx(pattern);
        if (!std::regex_match(cmp, rgx))
            throw parser_invalid_argument(detail::to_string("Value ", cmp, " did not match the pattern ", pattern, "."));
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must match the pattern '", pattern, "'.");
    }

private:
    //!\brief The pattern to match.
    std::string pattern;
};

/*!\brief A validator that checks if each value in a container satisfies a regex expression.
 * \ingroup argument_parser
 * \extends seqan3::regex_validator
 */
template <>
class regex_validator<std::vector<std::string>>
{
public:
    //!\brief Type of values that are tested by validator.
    using value_type = std::vector<std::string>;
    /*!\brief Constructing from a vector.
     * \param[in] pattern_ The vector of valid values to test.
     */
    regex_validator(std::string const & pattern_) :
        pattern{pattern_}
    {}

    /*!\brief Tests whether cmp lies inside values.
     * \param[in] cmp The value to validate.
     * \throws parser_invalid_argument
     */
    void operator()(value_type const & cmp) const
    {
        std::regex rgx(pattern);
        for (auto const & cmp_v : cmp)
            if (!std::regex_match(cmp_v, rgx))
                throw parser_invalid_argument(detail::to_string("Value ", cmp_v, " did not match the pattern ", pattern, "."));
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must match the pattern '", pattern, "'.");
    }

private:
    //!\brief The pattern to match.
    std::string pattern;
};

namespace detail
{

/*!\brief Validator that always returns true.
 * \ingroup argument_parser
 *
 * \details
 *
 * The default validator is needed to make the validator parameter of
 * argument_parser::add_option and argument_parser::add_option optional.
 */
template <typename option_value_type>
struct default_validator
{
    //!\brief Type of values that are tested by validator
    using value_type = option_value_type;

    //!\brief Value cmp always passes validation.
    bool operator()(option_value_type const & /*cmp*/) const noexcept
    {
        return true;
    }

    //!\brief Since no validation is happening the help message is empty.
    std::string get_help_page_message() const
    {
        return "";
    }
};

} // namespace detail

} // namespace seqan3
