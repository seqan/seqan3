// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides some standard validators for (positional) options.
 */

#pragma once

#include <regex>
#include <fstream>
#include <sstream>

#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/to_string.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/detail/safe_filesystem_entry.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\interface seqan3::validator <>
 * \brief The concept for option validators passed to add_option/positional_option.
 * \ingroup argument_parser
 *
 * \details
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
/*!\name Requirements for seqan3::validator
 * \brief You can expect these (meta-)functions on all types that implement seqan3::validator.
 * \{
 */
/*!\typedef     using option_value_type
 * \brief       The type of value on which the validator is called on.
 * \relates     seqan3::validator
 *
 * \details
 * \attention This is a concept requirement, not an actual typedef (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn              void operator()(option_value_type const & cmp) const
 * \brief           Validates the value 'cmp' and throws a seqan3::validation_error on failure.
 * \tparam          option_value_type The type of the value to be validated.
 * \param[in,out]   cmp The value to be validated.
 * \relates         seqan3::validator
 * \throws          seqan3::validation_error if value 'cmp' does not pass validation.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn              std::string get_help_page_message() const
 * \brief           Returns a message that can be appended to the (positional) options help page info.
 * \relates         seqan3::validator
 * \returns         A string that contains information on the performed validation.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
//!\cond
template <typename validator_type>
SEQAN3_CONCEPT validator = std::copyable<remove_cvref_t<validator_type>> &&
                           requires(validator_type validator,
                                    typename std::remove_reference_t<validator_type>::option_value_type value)
{
    typename std::remove_reference_t<validator_type>::option_value_type;

    { validator(value) } -> void;
    { validator.get_help_page_message() } -> std::string;
};
//!\endcond

/*!\brief A validator that checks whether a number is inside a given range.
 * \ingroup argument_parser
 * \implements seqan3::validator
 *
 * \details
 *
 * On construction, the validator must receive a maximum and a minimum number.
 * The class than acts as a functor, that throws a seqan3::validation_error
 * exception whenever a given value does not lie inside the given min/max range.
 *
 * \include test/snippet/argument_parser/validators_1.cpp
 */
class arithmetic_range_validator
{
public:
    //!\brief The type of value that this validator invoked upon.
    using option_value_type = double;

    /*!\brief The constructor.
     * \param[in] min_ Minimum set for the range to test.
     * \param[in] max_ Maximum set for the range to test.
     */
    arithmetic_range_validator(option_value_type const min_, option_value_type const max_) :
        min{min_}, max{max_}
    {}

    /*!\brief Tests whether cmp lies inside [`min`, `max`].
     * \param cmp The input value to check.
     * \throws seqan3::validation_error
     */
    void operator()(option_value_type const & cmp) const
    {
        if (!((cmp <= max) && (cmp >= min)))
            throw validation_error{detail::to_string("Value ", cmp, " is not in range [", min, ",", max, "].")};
    }

    /*!\brief Tests whether every element in \p range lies inside [`min`, `max`].
     * \tparam range_type The type of range to check; must model std::ranges::forward_range. The value type must model
     *                    seqan3::arithmetic.
     * \param  range      The input range to iterate over and check every element.
     * \throws seqan3::validation_error
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires arithmetic<std::ranges::range_value_t<range_type>>
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
    option_value_type min{};

    //!\brief Maximum of the range to test.
    option_value_type max{};
};

/*!\brief A validator that checks whether a value is inside a list of valid values.
 * \ingroup argument_parser
 * \implements seqan3::validator
 * \tparam option_value_t Must be a (container of) arithmetic type(s).
 *
 * \details
 *
 * On construction, the validator must receive a range or parameter pack of valid values.
 * The class than acts as a functor, that throws a seqan3::validation_error
 * exception whenever a given value is not in the given list.
 *
 * \note In order to simplify the chaining of validators, the option value type is deduced to `double` for ranges whose
 *       value type models seqan3::arithmetic, and to `std::string` if the ranges value type is convertible to it.
 *       Otherwise, the option value type is deduced to the value type of the range.
 *
 * \include test/snippet/argument_parser/validators_2.cpp
 */
template <typename option_value_t>
class value_list_validator
{
public:
    //!\brief Type of values that are tested by validator
    using option_value_type = option_value_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    value_list_validator() = default;                                         //!< Defaulted.
    value_list_validator(value_list_validator const &) = default;             //!< Defaulted.
    value_list_validator(value_list_validator &&) = default;                  //!< Defaulted.
    value_list_validator & operator=(value_list_validator const &) = default; //!< Defaulted.
    value_list_validator & operator=(value_list_validator &&) = default;      //!< Defaulted.
    ~value_list_validator() = default;                                        //!< Defaulted.

    /*!\brief Constructing from a range.
     * \tparam range_type The type of range; must model std::ranges::forward_range and value_list_validator::option_value_type
     *                    must be constructible from the rvalue reference type of the given range.
     * \param[in] rng The range of valid values to test.
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::constructible_from<option_value_type, std::ranges::range_rvalue_reference_t<range_type>>
    //!\endcond
    value_list_validator(range_type rng)
    {
        values.clear();
        std::ranges::move(std::move(rng), std::ranges::back_inserter(values));
    }

    /*!\brief Constructing from a parameter pack.
     * \tparam option_types The type of option values in the parameter pack; The value_list_validator::option_value_type must
     *                      be constructible from each type in the parameter pack.
     * \param[in] opts The parameter pack values.
     */
    template <typename ...option_types>
    //!\cond
        requires ((std::constructible_from<option_value_type, option_types> && ...))
    //!\endcond
    value_list_validator(option_types && ...opts)
    {
        (values.emplace_back(std::forward<option_types>(opts)), ...);
    }
    //!\}

    /*!\brief Tests whether cmp lies inside values.
     * \param cmp The input value to check.
     * \throws seqan3::validation_error
     */
    void operator()(option_value_type const & cmp) const
    {
        if (!(std::find(values.begin(), values.end(), cmp) != values.end()))
            throw validation_error{detail::to_string("Value ", cmp, " is not one of ", std::views::all(values), ".")};
    }

    /*!\brief Tests whether every element in \p range lies inside values.
     * \tparam range_type The type of range to check; must model std::ranges::forward_range.
     * \param  range      The input range to iterate over and check every element.
     * \throws seqan3::validation_error
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::convertible_to<std::ranges::range_value_t<range_type>, option_value_type>
    //!\endcond
    void operator()(range_type const & range) const
    {
        std::for_each(std::ranges::begin(range), std::ranges::end(range), [&] (auto cmp) { (*this)(cmp); });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("Value must be one of ", std::views::all(values), ".");
    }

private:

    //!\brief Minimum of the range to test.
    std::vector<option_value_type> values{};
};

/*!\brief Type deduction guides
 * \relates seqan3::value_list_validator
 * \{
 */
//!\brief Deduction guide for a parameter pack over an arithmetic type.
template <arithmetic ...option_types>
//!\cond
    requires !(detail::is_char_adaptation_v<option_types> || ...)
//!\endcond
value_list_validator(option_types...) -> value_list_validator<double>;

//!\brief Deduction guide for ranges over an arithmetic type.
template <typename range_type>
//!\cond
    requires std::ranges::forward_range<std::decay_t<range_type>> &&
             arithmetic<std::ranges::range_value_t<range_type>> &&
             (!seqan3::detail::is_char_adaptation_v<std::ranges::range_value_t<range_type>>)
//!\endcond
value_list_validator(range_type && rng) -> value_list_validator<double>;

//!\brief Given a parameter pack of types that are convertible to std::string, delegate to value type std::string.
template <typename ...option_types>
//!\cond
    requires (std::constructible_from<std::string, std::decay_t<option_types>> && ...)
//!\endcond
value_list_validator(option_types...) -> value_list_validator<std::string>;

//!\brief Deduction guide for ranges over a value type convertible to std::string.
template <typename range_type>
//!\cond
    requires (std::ranges::forward_range<std::decay_t<range_type>> &&
              std::constructible_from<std::string, std::ranges::range_value_t<range_type>>)
//!\endcond
value_list_validator(range_type && rng) -> value_list_validator<std::string>;

//!\brief Deduction guide for a parameter pack.
template <typename ...option_types>
value_list_validator(option_types...) -> value_list_validator<seqan3::pack_traits::front<option_types...>>;

//!\brief Deduction guide for ranges.
template <typename range_type>
//!\cond
    requires (std::ranges::forward_range<std::decay_t<range_type>>)
//!\endcond
value_list_validator(range_type && rng) -> value_list_validator<std::ranges::range_value_t<range_type>>;
//!\}

/*!\brief An abstract base class for the file and directory validators.
 * \ingroup argument_parser
 * \tparam file_t The type of the file to get the valid extensions for; `void` on default.
 *
 * \details
 *
 * This class provides a common interface for seqan3::input_file_validator and the seqan3::output_file_validator as
 * well as the seqan3::input_directory_validator and seqan3::output_directory_validator.
 *
 * The type can be further specialised for the seqan3::input_file_validator and the seqan3::output_file_validator
 * using the template argument to determine the valid extensions from the given file type.
 */
class file_validator_base
{
public:

    //!\brief Type of values that are tested by validator.
    using option_value_type = std::string;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    file_validator_base() = default;                                        //!< Defaulted.
    file_validator_base(file_validator_base const &) = default;             //!< Defaulted.
    file_validator_base(file_validator_base &&) = default;                  //!< Defaulted.
    file_validator_base & operator=(file_validator_base const &) = default; //!< Defaulted.
    file_validator_base & operator=(file_validator_base &&) = default;      //!< Defaulted.
    virtual ~file_validator_base() = default;                               //!< Virtual destructor.
    //!\}

    /*!\brief Tests if the given path is a valid input, respectively output, file or directory.
     * \param path The path to validate.
     *
     * \details
     *
     * This is a pure virtual function and must be overloaded by the derived class.
     */
    virtual void operator()(std::filesystem::path const & path) const = 0;

    /*!\brief Tests whether every path in list \p v passes validation. See operator()(option_value_type const & value)
     *        for further information.
     * \tparam range_type The type of range to check; must model std::ranges::forward_range and the value type must
     *                    be convertible to std::filesystem::path.
     * \param  v          The input range to iterate over and check every element.
     * \throws seqan3::validation_error
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::convertible_to<std::ranges::range_value_t<range_type>, std::filesystem::path const &>
    //!\endcond
    void operator()(range_type const & v) const
    {
         std::for_each(v.begin(), v.end(), [&] (auto cmp) { this->operator()(cmp); });
    }

protected:
    /*!\brief Validates the given filename path based on the specified extensions.
     * \param path The filename path.
     * \throws seqan3::validation_error if the specified extensions don't match the given path, or
     *         std::filesystem::filesystem_error on underlying OS API errors.
     */
    void validate_filename(std::filesystem::path const & path) const
    {
        // If no valid extensions are given we can safely return here.
        if (extensions.empty())
            return;

        // Check if extension is available.
        if (!path.has_extension())
            throw validation_error{detail::to_string("The given filename ", path.string(),
                                                            " has no extension. Expected one of the following valid"
                                                            " extensions:", extensions, "!")};

        // Drop the dot.
        std::string drop_less_ext = path.extension().string().substr(1);

        // Compares the extensions in lower case.
        auto case_insensitive_equal_to = [&] (std::string const & ext)
        {
            return std::ranges::equal(ext, drop_less_ext, [] (char const chr1, char const chr2)
                   {
                       return std::tolower(chr1) == std::tolower(chr2);
                   });
        };

        // Check if requested extension is present.
        if (std::ranges::find_if(extensions, case_insensitive_equal_to) == extensions.end())
        {
            throw validation_error{detail::to_string("Expected one of the following valid extensions: ",
                                                             extensions, "! Got ", drop_less_ext, " instead!")};
        }
    }

    /*!\brief Checks if the given path is readable.
     * \param path The path to check.
     * \returns `true` if readable, otherwise `false`.
     * \throws seqan3::validation_error if the path is not readable, or
     *         std::filesystem::filesystem_error on underlying OS API errors.
     */
    void validate_readability(std::filesystem::path const & path) const
    {
        // Check if input directory is readable.
        if (std::filesystem::is_directory(path))
        {
            std::error_code ec{};
            std::filesystem::directory_iterator{path, ec};  // if directory iterator cannot be created, ec will be set.
            if (static_cast<bool>(ec))
                throw validation_error{detail::to_string("Cannot read the directory ", path ,"!")};
        }
        else
        {
            // Must be a regular file.
            if (!std::filesystem::is_regular_file(path))
                throw validation_error{detail::to_string("Expected a regular file ", path, "!")};

            std::ifstream file{path};
            if (!file.is_open() || !file.good())
                throw validation_error{detail::to_string("Cannot read the file ", path, "!")};
        }
    }

    /*!\brief Checks if the given path is writable.
     * \param path The path to check.
     * \returns `true` if writable, otherwise `false`.
     * \throws seqan3::validation_error if the file could not be opened for writing, or
     *         std::filesystem::filesystem_error on underlying OS API errors.
     */
    void validate_writeability(std::filesystem::path const & path) const
    {
        std::ofstream file{path};
        detail::safe_filesystem_entry file_guard{path};

        bool is_open = file.is_open();
        bool is_good = file.good();
        file.close();

        if (!is_good || !is_open)
            throw validation_error{detail::to_string("Cannot write ", path, "!")};

        file_guard.remove();
    }

    //!\brief Returns the information of valid file extensions.
    std::string valid_extensions_help_page_message() const
    {
        if (extensions.empty())
            return "";
        else
            return detail::to_string(" Valid file extensions are: [", extensions | views::join(std::string{", "}), "].");
    }

    //!\brief Stores the extensions.
    std::vector<std::string> extensions{};
};

/*!\brief A validator that checks if a given path is a valid input file.
 * \ingroup argument_parser
 * \implements seqan3::validator
 * \tparam file_t The type of the file to get the valid extensions for; `void` on default.
 *
 * \details
 *
 * On construction, the validator can receive a list (std::vector over std::string) of valid file extensions.
 * The class acts as a functor that throws a seqan3::validation_error exception whenever a given filename's
 * extension (std::filesystem::path) is not in the given list of valid file extensions, if the file does not exist, or
 * if the file does not have the proper read permissions.
 *
 * \include test/snippet/argument_parser/validators_input_file.cpp
 *
 * The valid extensions can also be obtained from a seqan3 formatted file type, e.g. seqan3::sequence_input_file, if
 * it is given as template argument to this class. The following snippet demonstrates the different ways to instantiate
 * the seqan3::output_file_validator.
 *
 * \include test/snippet/argument_parser/validators_input_file_ext_from_file.cpp
 *
 * \note The validator works on every type that can be implicitly converted to std::filesystem::path.
 */
template <typename file_t = void>
class input_file_validator : public file_validator_base
{
public:

    static_assert(std::same_as<file_t, void> || detail::has_type_valid_formats<file_t>,
                  "Expected either a template type with a static member called valid_formats (a file type) or void.");

    // Import from base class.
    using typename file_validator_base::option_value_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */

    /*!\brief Default constructor.
     *
     * \details
     *
     * If the class' template argument `file_t` names a valid seqan3 file type that contains a
     * static member `valid_formats`, e.g. seqan3::sequence_input_file::valid_formats, then it generates the
     * list of valid extensions from this file. Otherwise the extensions list is empty.
     */
    input_file_validator()
    {
        if constexpr (!std::same_as<file_t, void>)
            file_validator_base::extensions = detail::valid_file_extensions<typename file_t::valid_formats>();
    }

    input_file_validator(input_file_validator const &) = default;             //!< Defaulted.
    input_file_validator(input_file_validator &&) = default;                  //!< Defaulted.
    input_file_validator & operator=(input_file_validator const &) = default; //!< Defaulted.
    input_file_validator & operator=(input_file_validator &&) = default;      //!< Defaulted.
    virtual ~input_file_validator() = default;                                //!< Virtual destructor.

    /*!\brief Constructs from a given collection of valid extensions.
     * \param[in] extensions The valid extensions to validate for.
     *
     * \details
     *
     * This constructor is only available if `file_t` does not name a valid seqan3 file type that contains a
     * static member `valid_formats`.
     */
    explicit input_file_validator(std::vector<std::string> extensions)
    //!\cond
        requires std::same_as<file_t, void>
    //!\endcond
        : file_validator_base{}
    {
        file_validator_base::extensions = std::move(extensions);
    }

    // Import base class constructor.
    using file_validator_base::file_validator_base;
    //!\}

    // Import the base::operator()
    using file_validator_base::operator();

    /*!\brief Tests whether path is an existing regular file and is readable.
     * \param file The input value to check.
     * \throws seqan3::validation_error if the validation process failed. Might be nested with
     *         std::filesystem::filesystem_error on unhandled OS API errors.
     */
    virtual void operator()(std::filesystem::path const & file) const override
    {
        try
        {
            if (!std::filesystem::exists(file))
                throw validation_error{detail::to_string("The file ", file, " does not exist!")};

            // Check if file is regular and can be opened for reading.
            validate_readability(file);

            // Check extension.
            validate_filename(file);
        }
        catch (std::filesystem::filesystem_error & ex)
        {
            std::throw_with_nested(validation_error{"Unhandled filesystem error!"});
        }
        catch (...)
        {
            std::rethrow_exception(std::current_exception());
        }
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return "The input file must exist and read permissions must be granted." +
               valid_extensions_help_page_message();
    }
};

/*!\brief A validator that checks if a given path is a valid output file.
 * \ingroup argument_parser
 * \implements seqan3::validator
 * \tparam file_t The type of the file to get the valid extensions for; `void` on default.
 *
 * \details
 *
 * On construction, the validator can receive a list (std::vector over std::string) of valid file extensions.
 * The class acts as a functor that throws a seqan3::validation_error exception whenever a given filename's
 * extension (sts::string) is not in the given list of valid file extensions, if the file already exist, or if the
 * parent path does not have the proper writer permissions.
 *
 * \include test/snippet/argument_parser/validators_output_file.cpp
 *
 * The valid extensions can also be obtained from a seqan3 formatted file type, e.g. seqan3::sequence_input_file, if
 * it is given as template argument to this class. The following snippet demonstrates the different ways to instantiate
 * the seqan3::output_file_validator.
 *
 * \include test/snippet/argument_parser/validators_output_file_ext_from_file.cpp
 *
 * \note The validator works on every type that can be implicitly converted to std::filesystem::path.
 */
template <typename file_t = void>
class output_file_validator : public file_validator_base
{
public:

    static_assert(std::same_as<file_t, void> || detail::has_type_valid_formats<file_t>,
                "Expected either a template type with a static member called valid_formats (a file type) or void.");

    // Import from base class.
    using typename file_validator_base::option_value_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */

    //!\copydoc seqan3::input_file_validator::input_file_validator()
    output_file_validator()
    {
        if constexpr (!std::same_as<file_t, void>)
            file_validator_base::extensions = detail::valid_file_extensions<typename file_t::valid_formats>();
    }

    output_file_validator(output_file_validator const &) = default;               //!< Defaulted.
    output_file_validator(output_file_validator &&) = default;                    //!< Defaulted.
    output_file_validator & operator=(output_file_validator const &) = default;   //!< Defaulted.
    output_file_validator & operator=(output_file_validator &&) = default;        //!< Defaulted.
    virtual ~output_file_validator() = default;                                   //!< Virtual Destructor.

    //!\copydoc seqan3::input_file_validator::input_file_validator(std::vector<std::string>)
    explicit output_file_validator(std::vector<std::string> extensions)
    //!\cond
        requires std::same_as<file_t, void>
    //!\endcond
        : file_validator_base{}
    {
        file_validator_base::extensions = std::move(extensions);
    }

    // Import base constructor.
    using file_validator_base::file_validator_base;
    //!\}

    // Import the base::operator()
    using file_validator_base::operator();

    /*!\brief Tests whether path is does not already exists and is writable.
     * \param file The input value to check.
     * \throws seqan3::validation_error if the validation process failed. Might be nested with
     *         std::filesystem::filesystem_error on unhandled OS API errors.
     */
    virtual void operator()(std::filesystem::path const & file) const override
    {
        try
        {
            if (std::filesystem::exists(file))
                throw validation_error{detail::to_string("The file ", file, " already exists!")};

            // Check if file has any write permissions.
            validate_writeability(file);

            validate_filename(file);
        }
        catch (std::filesystem::filesystem_error & ex)
        {
            std::throw_with_nested(validation_error{"Unhandled filesystem error!"});
        }
        catch (...)
        {
            std::rethrow_exception(std::current_exception());
        }
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return "The output file must not exist already and write permissions must be granted." +
               valid_extensions_help_page_message();
    }
};

/*!\brief A validator that checks if a given path is a valid input directory.
 * \ingroup argument_parser
 * \implements seqan3::validator
 *
 * \details
 *
 * The class acts as a functor that throws a seqan3::validation_error exception whenever a given directory
 * (std::filesystem::path) does not exist, the specified path is not a directory, or if the directory is not
 * readable.
 *
 * \include test/snippet/argument_parser/validators_input_directory.cpp
 *
 * \note The validator works on every type that can be implicitly converted to std::filesystem::path.
 */
class input_directory_validator : public file_validator_base
{
public:
    // Import from base class.
    using typename file_validator_base::option_value_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    input_directory_validator() = default;                                              //!< Defaulted.
    input_directory_validator(input_directory_validator const &) = default;             //!< Defaulted.
    input_directory_validator(input_directory_validator &&) = default;                  //!< Defaulted.
    input_directory_validator & operator=(input_directory_validator const &) = default; //!< Defaulted.
    input_directory_validator & operator=(input_directory_validator &&) = default;      //!< Defaulted.
    virtual ~input_directory_validator() = default;                                     //!< Virtual Destructor.

    // Import base constructor.
    using file_validator_base::file_validator_base;
    //!\}

    // Import the base::operator()
    using file_validator_base::operator();

    /*!\brief Tests whether path is an existing directory and is readable.
     * \param dir The input value to check.
     * \throws seqan3::validation_error if the validation process failed. Might be nested with
     *         std::filesystem::filesystem_error on unhandled OS API errors.
     */
    virtual void operator()(std::filesystem::path const & dir) const override
    {
        try
        {
            if (!std::filesystem::exists(dir))
                throw validation_error{detail::to_string("The directory ", dir, " does not exists!")};

            if (!std::filesystem::is_directory(dir))
                throw validation_error{detail::to_string("The path ", dir, " is not a directory!")};

            // Check if directory has any read permissions.
            validate_readability(dir);
        }
        catch (std::filesystem::filesystem_error & ex)
        {
            std::throw_with_nested(validation_error{"Unhandled filesystem error!"});
        }
        catch (...)
        {
            std::rethrow_exception(std::current_exception());
        }
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("An existing, readable path for the input directory.");
    }
};

/*!\brief A validator that checks if a given path is a valid output directory.
 * \ingroup argument_parser
 * \implements seqan3::validator
 *
 * \details
 *
 * The class acts as a functor that throws a seqan3::validation_error exception whenever a given path
 * (std::filesystem::path) is not writable. This can happen if either the parent path does not exists, or the
 * path doesn't have the proper write permissions.
 *
 * \include test/snippet/argument_parser/validators_output_directory.cpp
 *
 * \note The validator works on every type that can be implicitly converted to std::filesystem::path.
 */
class output_directory_validator : public file_validator_base
{
public:
    // Imported from base class.
    using typename file_validator_base::option_value_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    output_directory_validator() = default;                                               //!< Defaulted.
    output_directory_validator(output_directory_validator const &) = default;             //!< Defaulted.
    output_directory_validator(output_directory_validator &&) = default;                  //!< Defaulted.
    output_directory_validator & operator=(output_directory_validator const &) = default; //!< Defaulted.
    output_directory_validator & operator=(output_directory_validator &&) = default;      //!< Defaulted.
    virtual ~output_directory_validator() = default;                                      //!< Virtual Destructor.

    // Import base constructor.
    using file_validator_base::file_validator_base;
    //!\}

    // Import the base::operator().
    using file_validator_base::operator();

    /*!\brief Tests whether path is writable.
     * \param dir The input value to check.
     * \throws seqan3::validation_error if the validation process failed. Might be nested with
     *         std::filesystem::filesystem_error on unhandled OS API errors.
     */
    virtual void operator()(std::filesystem::path const & dir) const override
    {
        bool dir_exists = std::filesystem::exists(dir);
        // Make sure the created dir is deleted after we are done.
        std::error_code ec;
        std::filesystem::create_directory(dir, ec); // does nothing and is not treated as error if path already exists.
        // if error code was set or if dummy.txt could not be created within the output dir, throw an error.
        if (static_cast<bool>(ec))
            throw validation_error{detail::to_string("Cannot create directory: ", dir, "!")};

        try
        {
            if (!dir_exists)
            {
                detail::safe_filesystem_entry dir_guard{dir};
                validate_writeability(dir / "dummy.txt");
                dir_guard.remove_all();
            }
            else
            {
                validate_writeability(dir / "dummy.txt");
            }
        }
        catch (std::filesystem::filesystem_error & ex)
        {
            std::throw_with_nested(validation_error{"Unhandled filesystem error!"});
        }
        catch (...)
        {
            std::rethrow_exception(std::current_exception());
        }
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        return detail::to_string("A valid path for the output directory.");
    }
};

/*!\brief A validator that checks if a matches a regular expression pattern.
 * \ingroup argument_parser
 * \implements seqan3::validator
 *
 * \details
 *
 * On construction, the validator must receive a pattern for a regular expression.
 * The pattern variable will be used for constructing an std::regex and the
 * validator will call std::regex_match on the command line argument.
 * Note: A regex_match will only return true if the strings matches the pattern
 * completely (in contrast to regex_search which also matches substrings).
 *
 * The class than acts as a functor, that throws a seqan3::validation_error
 * exception whenever string does not match the pattern.
 *
 * \include test/snippet/argument_parser/validators_4.cpp
 */
class regex_validator
{
public:
    //!\brief Type of values that are tested by validator.
    using option_value_type = std::string;

    /*!\brief Constructing from a vector.
     * \param[in] pattern_ The pattern to match.
     */
    regex_validator(std::string const & pattern_) :
        pattern{pattern_}
    {}

    /*!\brief Tests whether cmp lies inside values.
     * \param[in] cmp The value to validate.
     * \throws seqan3::validation_error
     */
    void operator()(option_value_type const & cmp) const
    {
        std::regex rgx(pattern);
        if (!std::regex_match(cmp, rgx))
            throw validation_error{detail::to_string("Value ", cmp, " did not match the pattern ", pattern, ".")};
    }

    /*!\brief Tests whether every filename in list v matches the pattern.
     * \tparam range_type The type of range to check; must model std::ranges::forward_range and the value type must
     *                    be convertible to std::string.
     * \param  v          The input range to iterate over and check every element.
     * \throws seqan3::validation_error
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::convertible_to<std::ranges::range_value_t<range_type>, option_value_type const &>
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
 * \implements seqan3::validator
 * \tparam option_value_t Must be a (container of) arithmetic type(s).
 *
 * \details
 *
 * The default validator is needed to make the validator parameter of
 * argument_parser::add_option and argument_parser::add_option optional.
 */
template <typename option_value_t>
struct default_validator
{
    //!\brief Type of values that are tested by validator
    using option_value_type = option_value_t;

    //!\brief Value cmp always passes validation because the operator never throws.
    void operator()(option_value_t const & /*cmp*/) const noexcept
    {}

    //!\brief Since no validation is happening the help message is empty.
    std::string get_help_page_message() const
    {
        return "";
    }
};

/*!\brief A helper struct to chain validators recursively via the pipe operator.
 *\ingroup argument_parser
 *\implements seqan3::validator
 *
 *\details
 *
 * Note that both validators must operate on the same option_value_type in order to
 * avoid unexpected behaviour and ensure that the seqan3::argument_parser::add_option
 * call is well-formed. (add_option(val, ...., validator) requires
 * that val is of same type as validator::option_value_type).
 */
template <validator validator1_type, validator validator2_type>
//!\cond
    requires std::same_as<typename validator1_type::option_value_type, typename validator2_type::option_value_type>
//!\endcond
class validator_chain_adaptor
{
public:
    //!\brief The underlying type in both validators.
    using option_value_type = typename validator1_type::option_value_type;

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
        requires std::invocable<validator1_type, cmp_type const> && std::invocable<validator2_type, cmp_type const>
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
 * \ingroup argument_parser
 * \tparam validator1_type The type of the fist validator;
 *                         Must satisfy the seqan3::validator and the
 *                         same option_value_type as the second validator type.
 * \tparam validator2_type The type of the second validator;
 *                         Must satisfy the seqan3::validator and the
 *                         same option_value_type as the fist validator type.
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
 * seqan3::input_file_validator like this:
 *
 * \include test/snippet/argument_parser/validators_chaining.cpp
 *
 * You can chain as many validators as you want which will be evaluated one after
 * the other from left to right (first to last).
 */
template <validator validator1_type, validator validator2_type>
//!\cond
    requires std::same_as<typename std::remove_reference_t<validator1_type>::option_value_type,
                       typename std::remove_reference_t<validator2_type>::option_value_type>
//!\endcond
auto operator|(validator1_type && vali1, validator2_type && vali2)
{
    return detail::validator_chain_adaptor{std::forward<validator1_type>(vali1),
                                           std::forward<validator2_type>(vali2)};
}

} // namespace seqan3
