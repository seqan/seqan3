// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::default_printer.
 */

#pragma once

#include <iosfwd>
#include <tuple>
#include <utility>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

// clang-format off
/*!
 * \brief A tag that indicates that no printer was found for the given type.
 * \ingroup core_debug_stream
 */
struct no_printer_found{};
template <typename> struct advanceable_alignment_coordinate_printer {};
template <typename> struct alignment_matrix_printer {};
template <typename> struct alignment_printer {};
template <typename> struct alignment_result_printer {};
template <typename> struct alphabet_printer {};
template <typename> struct debug_stream_printer {};
template <typename> struct input_range_printer {};
template <typename> struct integer_sequence_printer {};
template <typename> struct integral_printer {};
template <typename> struct mask_printer {};
template <typename> struct optional_printer {};
template <typename> struct sequence_printer {};
template <typename> struct std_printer {};
template <typename> struct char_sequence_printer {};
template <typename> struct trace_directions_printer {};
template <typename> struct tuple_printer {};
// clang-format on

/*!
 * \interface seqan3::printable_with <>
 * \brief The concept for a printable object.
 *
 * A printable object is a printer that can print an argument to a stream.
 * The printer must be invocable with a stream and an argument.
 *
 * \tparam printer_t The type of the printer.
 * \tparam stream_t The type of the stream.
 * \tparam arg_t The type of the argument.
 * \ingroup core_debug_stream
 */
template <typename printer_t, typename stream_t, typename arg_t>
concept printable_with = std::invocable<printer_t, stream_t, arg_t>;

/*!
 * \brief The printer for standard output streams.
 *
 * The std_printer is used as a generic fallback to print regular types that can be printed to a
 * standard output stream, e.g. std::cout.
 *
 * \tparam type_t The type of the printable argument.
 * \ingroup core_debug_stream
 */
template <typename type_t>
    requires requires (std::ostream & cout, type_t const & value) {
                 {
                     cout << value
                 };
             }
struct std_printer<type_t>
{
    /*!
     * \brief The function call operator that prints the value to the stream.
     *
     * This function call operator prints the value to underlying stream of the seqan3::debug_stream.
     * This overload is only provided if the stream is a seqan3::debug_stream.
     *
     * \tparam stream_t The type of the stream.
     * \param[in,out] stream The stream to print to.
     * \param[in] value The value to print.
     */
    template <typename stream_t>
        requires requires (stream_t & stream) { stream.stream; }
    constexpr void operator()(stream_t & stream, type_t const & value) const
    {
        *stream.stream << value;
    }
};

/*!
 8 \brief The printer for integral types.
 *
 * The integral_printer is used to print integral types to a stream.
 *
 * \tparam integral_t The type of the integral.
 * \ingroup core_debug_stream
 */
template <typename integral_t>
    requires std::integral<std::remove_cvref_t<integral_t>>
struct integral_printer<integral_t>
{
    /*!
     * \brief The function call operator that prints the integral to the stream.
     * \tparam stream_t The type of the stream.
     * \param[in,out] stream The stream to print to.
     * \param[in] value The integral to print.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, integral_t const value) const
    {
        // note that we assume here that we always can print all std::integral's,
        // but this is not correct since std::cout << char32_t{5}; is not possible.
        // since char32_t is also an alphabet, we avoid infinite recursion here.
        if constexpr (printable_with<std_printer<integral_t>, stream_t &, integral_t>)
            std::invoke(std_printer<integral_t>{}, stream, value);
        else
            static_assert(std::same_as<integral_t, void>, "This type is not printable.");
    }
};

/*!
 * \brief The printer_order is a variadic template that defines the order of the printers.
 *
 * The printer_order is a variadic template that defines the order of the printers.
 * It is used to find the first valid printer among all given printers that can print the argument to the given stream.
 * The printer_order is used by the seqan3::default_printer.
 *
 * \tparam printer_templates_t The printer class templates that are used to print the arguments.
 * \ingroup core_debug_stream
 */
template <template <typename> typename... printer_templates_t>
struct printer_order
{
protected:
    /*!
     * \brief Find the index of the first printer that can print the argument.
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \tparam printers_t The printer class templates that are used to print the arguments.
     * \returns The index of the printer that can print the argument or `n`, where `n` is equal to the
     * number of printers, if no printer was found.
     */
    template <typename stream_t, typename arg_t, typename... printers_t>
    static constexpr std::ptrdiff_t find_index()
    {
        std::ptrdiff_t i = 0;
        return ((printable_with<printers_t, stream_t, arg_t> ? false : ++i) && ...) ? sizeof...(printer_templates_t) : i;
    }

    /*!
     * \brief The type of the printer that can print the argument to the stream.
     *
     * In case no valid printer is found for the argument and stream type the seqan3::no_printer_found is used.
     *
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \tparam i The index of the printer that can print the argument (default initialised to the index of the printer
     * in the list that can print the argument to the stream).
     */
    template <typename stream_t, typename arg_t,
              std::ptrdiff_t i = find_index<stream_t, arg_t, printer_templates_t<arg_t>...>()>
    using printer_for_t =
    // first: the index of the printer that can print the arguments or sizeof...(printers_t) if no printer was found
    // second: the tuple of instantiated printers extended with no_printer_found
        std::tuple_element_t<i, std::tuple<printer_templates_t<arg_t>..., no_printer_found>>;
};

/*!
 * \brief The default printer that is used by seqan3::debug_stream.
 *
 * This printer offers a function call operator if the type is printable by one of the printers defined
 * in the seqan3::printer_order.
 * If no valid printer is found the function call operator is not defined.
 *
 * \ingroup core_debug_stream
 */
struct default_printer :
    public printer_order<
            debug_stream_printer,
            alignment_result_printer,                 // type seqan3::alignment_result<>
            alignment_printer,                        // concept seqan3::tuple_like<>
            advanceable_alignment_coordinate_printer, // type seqan3::detail::advanceable_alignment_coordinate<>
            alignment_matrix_printer,                 // concept seqan3::detail::matrix<>
            trace_directions_printer,                 // type seqan3::detail::trace_directions
            mask_printer,                             // type seqan3::mask
            integral_printer,                         // concept std::integral
            // NOTE: alphabet_printer needs the integral_printer overload, otherwise it might have infinite recursion due to
            // char and uint being an alphabet
            alphabet_printer,         // concept seqan3::alphabet
            char_sequence_printer,    // concept std::range::input_range<> with char value_type
            integer_sequence_printer, // concept std::range::input_range<> with std::integral value_type
            sequence_printer, // concept seqan3::sequence<>, i.e. std::range::input_range<> with seqan3::alphabet value_type
            input_range_printer, // concept std::range::input_range<>
            optional_printer,    // type std::optional<> or std::nullopt_t
            tuple_printer,       // concept seqan3::tuple_like<>
            std_printer          // anything that can be printed by std::ostream
    >
{
public:

    /*!
     * \brief The function call operator that is only defined if the type is printable.
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \param[in,out] stream The stream to print to.
     * \param[in] arg The argument to print.
     */
    template <typename stream_t, typename arg_t>
        requires printable_with<printer_for_t<stream_t &, arg_t>, stream_t &, arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        using printer_t = printer_for_t<stream_t &, arg_t>;
        std::invoke(printer_t{}, stream, std::forward<arg_t>(arg));
    }
};

} // namespace seqan3
