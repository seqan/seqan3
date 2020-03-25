// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various utility functions.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <variant>

#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/filesystem>
#include <seqan3/std/iterator>

namespace seqan3::detail
{

//!\brief Base class to deduce the std::variant type from format tags.
//!\ingroup io
template <typename list_t, template <typename ...> typename output_t>
struct variant_from_tags;

//!\brief Transfers a list of format tags (`...ts`) onto a std::variant by specialising output_t with each.
//!\ingroup io
template <template <typename...> typename output_t, typename ...ts>
struct variant_from_tags<type_list<ts...>, output_t>
{
    //!\brief The type of std::variant.
    using type = std::variant<output_t<ts>...>;
};

/*!\brief Write `"\n"` or `"\r\n"` to the stream iterator, depending on arguments.
 * \tparam  it_t Type of the iterator; must satisfy std::output_Iterator with `char`.
 * \param     it The iterator.
 * \param add_cr Whether to add carriage return, too.
 * \ingroup io
 */
template <std::output_iterator<char> it_t>
constexpr void write_eol(it_t & it, bool const add_cr)
{
    if (add_cr)
        it = '\r';

    it = '\n';
}

/*!\brief Sets the file format according to the file name extension.
 * \ingroup io
 * \tparam     format_variant_type The variant type of the format to set.
 * \param[out] format              The format to set.
 * \param[in]  file_name           The file name to extract the extension from.
 *
 * \throws seqan3::unhandled_extension_error If the extension in file_name does
 *         not occur in any valid extensions of the formats specified in the
 *         \p format_variant_type template argument list.
 */
template <typename format_variant_type>
void set_format(format_variant_type & format,
                std::filesystem::path const & file_name)
{
    using valid_formats = detail::transfer_template_args_onto_t<format_variant_type, type_list>;

    bool format_found = false;
    std::string extension = file_name.extension().string();
    if (extension.size() > 1)
    {
        extension = extension.substr(1); // drop leading "."
        detail::for_each<valid_formats>([&] (auto fmt)
        {
            using fm_type = typename decltype(fmt)::type; // remove type_identity wrapper

            for (auto const & ext : fm_type::file_extensions)
            {
                if (std::ranges::equal(ext, extension))
                {
                    format = fm_type{};
                    format_found = true;
                    return;
                }
            }
        });
    }

    if (!format_found)
        throw unhandled_extension_error("No valid format found for this extension.");
}

/*!\brief Helper function to determine if all types in a format type list have a static member `file_extensions`.
 * \tparam list_t The type of the template parameter list.
 * \returns `true` if `type::file_extensions` for all expanded types of `list_t` is valid, otherwise `false`.
 */
template <typename list_t>
inline constexpr bool has_member_file_extensions = false;

//!\cond
template <template <typename ...> typename list_t, typename ...ts>
    requires (requires { ts::file_extensions; }, ..., true)
inline constexpr bool has_member_file_extensions<list_t<ts...>> = true;
//!\endcond

/*!\brief Helper function to determine if a type has a static member `valid_formats`.
 * \tparam query_t The type to query.
 * \returns `true` if `query_t::valid_formats` is valid, otherwise `false`.
 */
template <typename query_t>
inline constexpr bool has_type_valid_formats = false;

//!\cond
template <typename query_t>
    requires requires { typename query_t::valid_formats; }
inline constexpr bool has_type_valid_formats<query_t> = true;
//!\endcond

/*!\brief Returns a list of valid file extensions.
 * \ingroup io
 * \tparam formats_t The list of formats to parse; seqan3::detail::all_formats_have_file_extensions must return `true`.
 * \returns `std::vector<std::string>` with all valid file extensions specified by `valid_formats`.
 *
 * \details
 *
 * ### Complexity
 *
 * Linear in the number of file extensions.
 *
 * ### Thread-safety
 *
 * Thread-safe.
 *
 * ### Exception
 *
 * Strong exception guarantee. No input is modified. Might throw std::bad_alloc.
 */
template <typename formats_t>
inline std::vector<std::string> valid_file_extensions()
{
    static_assert(has_member_file_extensions<formats_t>,
                  "Expects that all formats have a static member file_extensions storing the extensions in a range");

    std::vector<std::string> extensions;
    detail::for_each<formats_t>([&extensions] (auto t_identity)
    {
        using format_t = typename decltype(t_identity)::type;
        std::ranges::copy(format_t::file_extensions, std::ranges::back_inserter(extensions));
    });

    return extensions;
}
} // namespace seqan3::detail
