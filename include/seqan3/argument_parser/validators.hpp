// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides some standard validators for (positional) options.
 */

#pragma once

#include <seqan3/io/detail/misc.hpp>

#include <sharg/validators.hpp>

namespace seqan3
{

using namespace sharg;

/*!\brief A validator that checks if a given path is a valid input file.
 * \ingroup argument_parser
 * \implements seqan3::validator
 * \tparam file_t The type of the file to get the valid extensions for; `void` on default.
 *
 * \details
 *
 * The valid extensions can also be obtained from a seqan3 formatted file type, e.g. seqan3::sequence_input_file, if
 * it is given as template argument to this class. The following snippet demonstrates the different ways to instantiate
 * the seqan3::input_file_validator.
 *
 * \include test/snippet/argument_parser/validators_input_file_ext_from_file.cpp
 */
template <typename file_t = void>
class input_file_validator /* \cond */ : public sharg::input_file_validator /* \endcond */
{
public:
    static_assert(std::same_as<file_t, void> || detail::has_type_valid_formats<file_t>,
                  "Expected either a template type with a static member called valid_formats (a file type) or void.");

    // Import from base class.
    using typename sharg::input_file_validator::option_value_type;

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
    input_file_validator() : sharg::input_file_validator{default_extensions()}
    {}

    //!\cond
    // the documentation is included by the sharg parser.
    input_file_validator(input_file_validator const &) = default;
    input_file_validator(input_file_validator &&) = default;
    input_file_validator & operator=(input_file_validator const &) = default;
    input_file_validator & operator=(input_file_validator &&) = default;
    virtual ~input_file_validator() = default;
    //!\endcond

    /*!\brief Constructs from a given collection of valid extensions.
     * \param[in] extensions The valid extensions to validate for.
     *
     * \details
     *
     * This constructor is only available if `file_t` does not name a valid seqan3 file type that contains a
     * static member `valid_formats`.
     */
    explicit input_file_validator(std::vector<std::string> extensions)
        requires std::same_as<file_t, void>
    //!\endcond
    : sharg::input_file_validator{std::move(extensions)}
    {
        // sharg::input_file_validator::extensions = std::move(extensions);
    }

    // Import base class constructor.
    using sharg::input_file_validator::input_file_validator;
    //!\}

    /*!\brief The default extensions of `file_t`.
     * \returns A list of default extensions for `file_t`, will be empty if `file_t` is `void`.
     *
     * \details
     *
     * If `file_t` does name a valid seqan3 file type that contains a static member `valid_formats` returns the
     * extensions of that `file_t` type. Otherwise returns an empty list.
     */
    static std::vector<std::string> default_extensions()
    {
        if constexpr (!std::same_as<file_t, void>)
            return detail::valid_file_extensions<typename file_t::valid_formats>();
        return {};
    }

    // Import the base::operator()
    using sharg::input_file_validator::operator();
};

//!\cond
// guides are only needed because we inherit from the sharg::input_validator that does not have a template
template <typename file_t = void>
input_file_validator() -> input_file_validator<file_t>;
//!\endcond

/*!\brief A validator that checks if a given path is a valid output file.
 * \ingroup argument_parser
 * \implements seqan3::validator
 * \tparam file_t The type of the file to get the valid extensions for; `void` on default.
 *
 * \details
 *
 * The valid extensions can also be obtained from a seqan3 formatted file type, e.g. seqan3::sequence_input_file, if
 * it is given as template argument to this class. The following snippet demonstrates the different ways to instantiate
 * the seqan3::output_file_validator.
 *
 * \include test/snippet/argument_parser/validators_output_file_ext_from_file.cpp
 */
template <typename file_t = void>
class output_file_validator /* \cond */ : public sharg::output_file_validator /* \endcond */
{
public:
    static_assert(std::same_as<file_t, void> || detail::has_type_valid_formats<file_t>,
                  "Expected either a template type with a static member called valid_formats (a file type) or void.");

    // Import from base class.
    using typename sharg::output_file_validator::option_value_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */

    //!\copydoc seqan3::input_file_validator::input_file_validator()
    output_file_validator() :
        sharg::output_file_validator{output_file_open_options::create_new, this->default_extensions()}
    {}

    //!\cond
    // the documentation is included by the sharg parser.
    output_file_validator(output_file_validator const &) = default;             //!< Defaulted.
    output_file_validator(output_file_validator &&) = default;                  //!< Defaulted.
    output_file_validator & operator=(output_file_validator const &) = default; //!< Defaulted.
    output_file_validator & operator=(output_file_validator &&) = default;      //!< Defaulted.
    virtual ~output_file_validator() = default;                                 //!< Virtual destructor.
    //!\endcond

    /*!\brief Constructs from a given overwrite mode and a list of valid extensions.
     * \param[in] mode A seqan3::output_file_open_options indicating whether the validator throws if a file already
     *                 exists.
     */
    explicit output_file_validator(output_file_open_options const mode) :
        sharg::output_file_validator{mode, this->default_extensions()}
    {}

    // Import base constructor.
    using sharg::output_file_validator::output_file_validator;
    //!\}

    /*!\brief The default extensions of `file_t`.
     * \returns A list of default extensions for `file_t`, will be empty if `file_t` is `void`.
     *
     * \details
     *
     * If `file_t` does name a valid seqan3 file type that contains a static member `valid_formats` returns the
     * extensions of that `file_t` type. Otherwise returns an empty list.
     */
    static std::vector<std::string> default_extensions()
    {
        if constexpr (!std::same_as<file_t, void>)
            return detail::valid_file_extensions<typename file_t::valid_formats>();
        return {};
    }

    // Import the base::operator()
    using sharg::output_file_validator::operator();
};

//!\cond
// guides are only needed because we inherit from the sharg::input_validator that does not have a template
template <typename file_t = void>
output_file_validator() -> output_file_validator<file_t>;

template <typename file_t = void>
output_file_validator(output_file_open_options) -> output_file_validator<file_t>;

template <typename file_t = void>
output_file_validator(output_file_open_options, std::vector<std::string>) -> output_file_validator<file_t>;
//!\endcond

} // namespace seqan3
