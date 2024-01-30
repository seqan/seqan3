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

template <typename printer_t>
concept printable = requires () {
                        {
                            printer_t::print
                        };
                    };

template <typename type_t>
    requires requires (std::ostream & cout, type_t value) {
                 {
                     cout << value
                 };
             }
struct std_printer<type_t>
{
    static constexpr auto print = [](auto & s, auto && value)
    {
        *s.stream << std::forward<decltype(value)>(value);
    };
};

template <typename integral_t>
    requires std::integral<std::remove_cvref_t<integral_t>>
struct integral_printer<integral_t>
{
    static constexpr auto print = [](auto & s, auto const value)
    {
        // note that we assume here that we always can print all std::integral's,
        // but this is not correct since std::cout << char32_t{5}; is not possible.
        // since char32_t is also an alphabet, we avoid infinite recursion here.
        if constexpr (printable<std_printer<integral_t>>)
            std_printer<integral_t>::print(s, value);
        else
            static_assert(std::same_as<integral_t, void>, "This type is not printable.");
    };
};

template <template <typename> typename... printers_t>
struct printer_order
{
    template <typename type_t, typename... types>
    static constexpr std::ptrdiff_t find_index()
    {
        std::ptrdiff_t i = 0;
        return ((printable<types> ? false : ++i) && ...) ? -1 : i;
    }

    template <typename type_t, std::ptrdiff_t i = find_index<type_t, printers_t<type_t>...>()>
    using printer_for_t =
        std::tuple_element_t<i == -1 ? sizeof...(printers_t) : i, std::tuple<printers_t<type_t>..., no_printer_found>>;

    // std::tuple<...> list of all printers that can print type_t
    // NOTE: printer_order<...> might be nicer but is more complicated to write
    template <typename type_t>
    using printers_for_t = decltype(std::tuple_cat(
        std::conditional_t<printable<printers_t<type_t>>, std::tuple<printers_t<type_t>>, std::tuple<>>{}...));

    template <typename type_t>
    static constexpr bool is_printable = printable<printer_for_t<type_t>>;

    static constexpr auto print = [](auto & s, auto && arg)
    {
        using printer_t = printer_for_t<decltype(arg)>;
        printer_t::print(s, std::forward<decltype(arg)>(arg));
    };
};

struct default_printer :
    printer_order<
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
{};

} // namespace seqan3
