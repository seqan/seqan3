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

template <typename printer_t, typename stream_t, typename arg_t>
concept printable_with = std::invocable<printer_t, stream_t, arg_t>;

template <typename type_t>
    requires requires (std::ostream & cout, type_t const & value) {
                 {
                     cout << value
                 };
             }
struct std_printer<type_t>
{
    static constexpr auto print = [](auto & s, auto && value)
    {
        std_printer<type_t> printer{};
        std::invoke(printer, s, value);
    };

    template <typename stream_type>
    constexpr void operator()(stream_type & stream, type_t const & value) const
    {
        *stream.stream << value;
    }
};

template <typename integral_t>
    requires std::integral<std::remove_cvref_t<integral_t>>
struct integral_printer<integral_t>
{
    static constexpr auto print = [](auto & s, auto const value)
    {
        integral_printer<std::remove_cvref_t<integral_t>> printer{};
        std::invoke(printer, s, value);
    };

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

template <template <typename> typename... printer_templates_t>
struct printer_order
{
protected:
    template <typename stream_t, typename arg_t, typename... printers_t>
    static constexpr std::ptrdiff_t find_index()
    {
        std::ptrdiff_t i = 0;
        return ((printable_with<printers_t, stream_t, arg_t> ? false : ++i) && ...) ? -1 : i;
    }

    template <typename stream_t, typename arg_t,
              std::ptrdiff_t i = find_index<stream_t, arg_t, printer_templates_t<arg_t>...>()>
    using printer_for_t =
    // first: the index of the printer that can print the arguments or sizeof...(printers_t) if no printer was found
    // second: the tuple of instantiated printers extended with no_printer_found
        std::tuple_element_t<i == -1 ? sizeof...(printer_templates_t) : i,
                             std::tuple<printer_templates_t<arg_t>..., no_printer_found>>;
};

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

    template <typename stream_t, typename arg_t>
        requires printable_with<printer_for_t<stream_t &, arg_t>, stream_t &, arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        using printer_t = printer_for_t<stream_t &, arg_t>;
        printer_t printer{};
        std::invoke(printer, stream, std::forward<arg_t>(arg));
    }
};

} // namespace seqan3
