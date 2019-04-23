// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Contains some standard validators for (positional) options.
 */

#pragma once

#include <regex>
#include <sstream>

#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/to_string.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/to_lower.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

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
/*!\name Requirements for seqan3::validator_concept
 * \brief You can expect these (meta-)functions on all types that implement seqan3::validator_concept.
 * \{
 */
/*!\typedef     using value_type
 * \brief       The type of value on which the validator is called on.
 * \relates     seqan3::validator_concept
 *
 * \details
 * \attention This is a concept requirement, not an actual typedef (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn              void operator()(value_type const & cmp) const
 * \brief           Validates the value 'cmp' and throws a seqan3::validation_error on failure.
 * \tparam          value_type The type of the value to be validated.
 * \param[in,out]   cmp The value to be validated.
 * \relates         seqan3::validator_concept
 * \throws          seqan3::validation_error if value 'cmp' does not pass validation.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn              std::string get_help_page_message() const
 * \brief           Returns a message that can be appended to the (positional) options help page info.
 * \relates         seqan3::validator_concept
 * \returns         A string that contains information on the performed validation.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
//!\cond
template <typename validator_type>
SEQAN3_CONCEPT validator_concept = std::Copyable<remove_cvref_t<validator_type>> &&
                                   requires(validator_type validator,
                                            typename std::remove_reference_t<validator_type>::value_type value)
{
    typename std::remove_reference_t<validator_type>::value_type;

    { validator(value) } -> void;
    { validator.get_help_page_message() } -> std::string;
};
//!\endcond

/*!\brief A validator that checks whether a number is inside a given range.
 * \ingroup argument_parser
 * \implements seqan3::validator_concept
 *
 * \tparam option_value_type Must be a (container of) arithmetic type(s).
 *
 * \details
 *
 * On construction, the validator must receive a maximum and a minimum number.
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever a given value does not lie inside the given min/max range.
 *
 * \snippet test/snippet/argument_parser/validators_1.cpp usage
 */
class arithmetic_range_validator
{
public:
    //!\brief The type of value that this validator invoked upon.
    using value_type = double;

    /*!\brief The constructor.
     * \param[in] min_ Minimum set for the range to test.
     * \param[in] max_ Maximum set for the range to test.
     */
    arithmetic_range_validator(value_type const min_, value_type const max_) :
        min{min_}, max{max_}
    {}

    /*!\brief Tests whether cmp lies inside [`min`, `max`].
     * \param cmp The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(value_type const & cmp) const
    {
        if (!((cmp <= max) && (cmp >= min)))
            throw parser_invalid_argument(detail::to_string("Value ", cmp, " is not in range [", min, ",", max, "]."));
    }

    /*!\brief Tests whether every element in \p range lies inside [`min`, `max`].
     * \tparam range_type The type of range to check; must model std::ranges::ForwardRange. The value type must model
     *                    seqan3::Arithmetic.
     * \param  range      The input range to iterate over and check every element.
     * \throws parser_invalid_argument
     */
    template <std::ranges::ForwardRange range_type>
    //!\cond
        requires Arithmetic<value_type_t<range_type>>
    //!\endcond
    void operator()(range_type const & range) const
    {
        std::for_each(range.begin(), range.end(), [&] (auto cmp) { (*this)(cmp); });
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

/*!\brief A validator that checks whether a value is inside a list of valid values.
 * \ingroup argument_parser
 * \implements seqan3::validator_concept
 *
 * \details
 *
 * On construction, the validator must receive a list (vector) of valid values.
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever a given value is not in the given list.
 *
 * \snippet test/snippet/argument_parser/validators_2.cpp usage
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
    value_list_validator(std::vector<value_type> const & v) :
        values{v}
    {}

    /*!\brief Constructing from an initializer_list.
     * \param[in] v The initializer_list of valid values to test.
     */
    value_list_validator(std::initializer_list<value_type> const & v) :
        values{v}
    {}

    /*!\brief Tests whether cmp lies inside values.
     * \param cmp The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(value_type const & cmp) const
    {
        if (!(std::find(values.begin(), values.end(), cmp) != values.end()))
            throw parser_invalid_argument(detail::to_string("Value ", cmp, " is not one of ", std::view::all(values), "."));
    }

    /*!\brief Tests whether every element in \p range lies inside values.
     * \tparam range_type The type of range to check; must model std::ranges::ForwardRange.
     * \param  range      The input range to iterate over and check every element.
     * \throws parser_invalid_argument
     */
    template <std::ranges::ForwardRange range_type>
    //!\cond
        requires std::ConvertibleTo<value_type_t<range_type>, value_type const &>
    //!\endcond
    void operator()(range_type const & range) const
    {
        std::for_each(range.begin(), range.end(), [&] (auto cmp) { (*this)(cmp); });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must be one of ", std::view::all(values), ".");
    }

private:

    //!\brief Minimum of the range to test.
    std::vector<value_type> values;
};

/*!\brief Type deduction guides
 * \relates seqan3::value_list_validator
 * \{
 */
template <Arithmetic option_value_type>
value_list_validator(std::vector<option_value_type> const & v) -> value_list_validator<double>;

template <Arithmetic option_value_type>
value_list_validator(std::initializer_list<option_value_type> const & v) -> value_list_validator<double>;

value_list_validator(std::vector<const char *> const & v) -> value_list_validator<std::string>;

value_list_validator(std::initializer_list<const char *> const & v) -> value_list_validator<std::string>;
//!\}

/*!\brief A validator that checks if a filenames has one of the valid extensions.
 * \ingroup argument_parser
 * \implements seqan3::validator_concept
 *
 * \details
 *
 * On construction, the validator must receive a list (vector) of valid file extensions.
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever a given filename (string) is not in the given extension list.
 *
 * \snippet test/snippet/argument_parser/validators_3.cpp usage
 *
 * \note The validator works on every type that can be implicitly cast to std::filesystem::path.
 */
class file_ext_validator
{
public:
    //!\brief Type of values that are tested by validator
    using value_type = std::string;

    /*!\brief Constructing from a vector.
     * \param[in] v The vector of valid file extensions to test (e.g. {"fa", "fasta"}).
     * \param[in] c Case sensitivity flag. Set true for case sensitivity. Default: false (case insensitive).
     *
     * For case insensitivity, everything is converted to lower case characters.
     */
    file_ext_validator(std::vector<std::string> const & v, bool const c = false) :
        case_sensitive{c}
    {
        extensions = c ? v : std::vector<std::string>{v | view::to_lower};
    }

    /*!\brief Tests whether the filepath \p path ends with a valid extension.
     * \param path The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(std::filesystem::path const & path) const
    {
        std::string ext{path.extension().string()};
        ext = ext.substr(std::min(1, std::max(0, static_cast<int>(ext.size()) - 1))); // drop '.' if ext is non-empty

        // extensions were transformed to lower case during construction,so do the same for input path
        if (!case_sensitive)
            ext = std::string{ext | view::to_lower};

        if (!(std::find(extensions.begin(), extensions.end(), ext) != extensions.end()))
            throw parser_invalid_argument(detail::to_string("Extension ", ext, " is not one of ",
                                                            std::view::all(extensions), "."));
    }

    /*!\brief Tests whether every value of v lies inside extensions.
     * \tparam range_type The type of range to check; must model std::ranges::ForwardRange and the value type must
     *                    be convertible to std::filesystem::path.
     * \param  v          The input range to iterate over and check every element.
     * \throws parser_invalid_argument
     */
    template <std::ranges::ForwardRange range_type>
    //!\cond
        requires std::ConvertibleTo<value_type_t<range_type>, std::filesystem::path const &>
    //!\endcond
    void operator()(range_type const & v) const
    {
        std::for_each(v.begin(), v.end(), [&] (auto cmp) { (*this)(cmp); });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("File name extension must be one of ", std::view::all(extensions), ".");
    }

private:
    //!\brief Stores valid file extensions.
    std::vector<std::string> extensions;

    //!\brief True if file extension is case sensitive
    bool case_sensitive;
};

/*!\brief A validator that checks if a path (file or directory) exists.
 * \ingroup argument_parser
 *
 * \details
 *
 * The struct then acts as a functor that throws a seqan3::parser_invalid_argument
 * exception whenever a given path (file or directory) does not exist.
 *
 * \snippet test/snippet/argument_parser/validators_path_existence.cpp usage
 */
class path_existence_validator
{
public:
    //!\brief Type of values that are tested by validator
    using value_type = std::string;

    /*!\brief Tests whether path (file or directory) exists.
     * \param path The input value to check.
     * \throws parser_invalid_argument
     */
    void operator()(std::filesystem::path const & path) const
    {
        if (!(std::filesystem::exists(path)))
            throw parser_invalid_argument(detail::to_string("The file or directory ", path, " does not exist."));
    }

    /*!\brief Tests whether every path (file or directory) in list \p v exists.
     * \tparam range_type The type of range to check; must model std::ranges::ForwardRange and the value type must
     *                    be convertible to std::filesystem::path.
     * \param  v          The input range to iterate over and check every element.
     * \throws parser_invalid_argument
     */
    template <std::ranges::ForwardRange range_type>
    //!\cond
        requires std::ConvertibleTo<value_type_t<range_type>, std::filesystem::path const &>
    //!\endcond
    void operator()(range_type const & v) const
    {
         std::for_each(v.begin(), v.end(), [&] (auto cmp) { (*this)(cmp); });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("The file or directory is checked for existence.");
    }
};

/*!\brief A validator that checks if a matches a regular expression pattern.
 * \ingroup argument_parser
 * \implements seqan3::validator_concept
 *
 * \details
 *
 * On construction, the validator must receive a pattern for a regular expression.
 * The pattern variable will be used for constructing an std::regex and the
 * validator will call std::regex_match on the command line argument.
 * Note: A regex_match will only return true if the strings matches the pattern
 * completely (in contrast to regex_search which also matches substrings).
 *
 * The struct than acts as a functor, that throws a seqan3::parser_invalid_argument
 * exception whenever string does not match the pattern.
 *
 * \snippet test/snippet/argument_parser/validators_4.cpp usage
 */
class regex_validator
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

    /*!\brief Tests whether every filename in list v matches the pattern.
     * \tparam range_type The type of range to check; must model std::ranges::ForwardRange and the value type must
     *                    be convertible to std::string.
     * \param  v          The input range to iterate over and check every element.
     * \throws parser_invalid_argument
     */
    template <std::ranges::ForwardRange range_type>
    //!\cond
        requires std::ConvertibleTo<value_type_t<range_type>, value_type const &>
    //!\endcond
    void operator()(range_type const & v) const
    {
         std::for_each(v.begin(), v.end(), [&] (auto cmp) { (*this)(cmp); });
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
 * \implements seqan3::validator_concept
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

    //!\brief Value cmp always passes validation because the operator never throws.
    void operator()(option_value_type const & /*cmp*/) const noexcept
    {}

    //!\brief Since no validation is happening the help message is empty.
    std::string get_help_page_message() const
    {
        return "";
    }
};

/*!\brief A helper struct to chain validators recursively via the pipe operator.
 *\ingroup argument_parser
 *\implements seqan3::validator_concept
 *
 *\details
 *
 * Note that both validators must operate on the same value_type in order to
 * avoid unexpected behaviour and ensure that the seqan3::argument_parser::add_option
 * call is well-formed. (add_option(val, ...., validator) requires
 * that val is of same type as validator::value_type).
 */
template <validator_concept validator1_type, validator_concept validator2_type>
//!\cond
    requires std::Same<typename validator1_type::value_type, typename validator2_type::value_type>
//!\endcond
class validator_chain_adaptor
{
public:
    //!\brief The underlying type in both validators.
    using value_type = typename validator1_type::value_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    validator_chain_adaptor() = delete;                                                //!< Deleted.
    validator_chain_adaptor(validator_chain_adaptor const & pf) = default;             //!< Defaulted.
    validator_chain_adaptor & operator=(validator_chain_adaptor const & pf) = default; //!< Defaulted.
    validator_chain_adaptor(validator_chain_adaptor &&) = default;                     //!< Defaulted.
    validator_chain_adaptor & operator=(validator_chain_adaptor &&) = default;         //!< Defaulted.

    /*!\brief Constructing from two validators.
     * \param[in] vali1_ Some validator to be chained to vali2_.
     * \param[in] vali2_ Another validator to be chained to vali1_.
     */
    validator_chain_adaptor(validator1_type vali1_, validator2_type vali2_) :
        vali1{std::move(vali1_)}, vali2{std::move(vali2_)}
    {}

    //!\brief The destructor.
    ~validator_chain_adaptor() = default;
    //!\}

    /*!\brief Calls the operator() of each validator on the value cmp.
     * \tparam cmp_type The type of value to validate; must be invokable with each of the validator members.
     * \param[in] cmp   The value to validate.
     *
     * This function delegates to the validation of both of the chained validators
     * by calling their operator() one after the other. The behaviour depends on
     * the chained validators which may throw on input error.
     */
    template <typename cmp_type>
    //!\cond
        requires std::Invocable<validator1_type, cmp_type const> && std::Invocable<validator2_type, cmp_type const>
    //!\endcond
    void operator()(cmp_type const & cmp) const
    {
        vali1(cmp);
        vali2(cmp);
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string(vali1.get_help_page_message(), " ", vali2.get_help_page_message());
    }

private:
    //!\brief The first validator in the chain.
    validator1_type vali1;
    //!\brief The second validator in the chain.
    validator2_type vali2;
};

} // namespace detail

/*!\brief Enables the chaining of validators.
 *!\ingroup argument_parser
 * \tparam validator1_type The type of the fist validator;
 *                         Must satisfy the seqan3::validator_concept and the
 *                         same value_type as the second validator type.
 * \tparam validator2_type The type of the second validator;
 *                         Must satisfy the seqan3::validator_concept and the
 *                         same value_type as the fist validator type.
 * \param[in] vali1 The first validator to chain.
 * \param[in] vali2 The second validator to chain.
 * \returns A new validator that tests a value for both vali1 and vali2.
 *
 * \details
 *
 * The pipe operator is the AND operation for two validators, which means that a
 * value must pass both validators in order to be accepted by the new validator.
 *
 * For example you may want a file name that only accepts absolute paths but
 * also must have one out of some given file extensions.
 * For this purpose you can chain a seqan3::regex_validator to a
 * seqan3::file_ext_validator like this:
 *
 * \include test/snippet/argument_parser/validators_chaining.cpp
 *
 * You can chain as many validators as you want which will be evaluated one after
 * the other from left to right (first to last).
 */
template <validator_concept validator1_type, validator_concept validator2_type>
//!\cond
    requires std::Same<typename std::remove_reference_t<validator1_type>::value_type,
                       typename std::remove_reference_t<validator2_type>::value_type>
//!\endcond
auto operator|(validator1_type && vali1, validator2_type && vali2)
{
    return detail::validator_chain_adaptor{std::forward<validator1_type>(vali1),
                                           std::forward<validator2_type>(vali2)};
}

} // namespace seqan3
