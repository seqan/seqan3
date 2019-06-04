// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various utility functions.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <variant>

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/filesystem>
#include <seqan3/std/iterator>

namespace seqan3::detail
{

//!\brief Base class to deduce the std::variant type from format tags.
//!\ingroup io
template <typename list_t, template <typename...> typename output_t>
struct variant_from_tags;

//!\brief Transfers a list of format tags (`...ts`) onto a std::variant by specialising output_t with each.
//!\ingroup io
template <template <typename...> typename output_t, typename ...ts>
struct variant_from_tags<meta::list<ts...>, output_t>
{
    //!\brief the type of std::variant.
    using type = std::variant<output_t<ts>...>;
};

/*!\brief Write `"\n"` or `"\r\n"` to the stream iterator, depending on arguments.
 * \tparam  it_t Type of the iterator; must satisfy std::output_Iterator with `char`.
 * \param     it The iterator.
 * \param add_cr Whether to add carriage return, too.
 * \ingroup io
 */
template <std::OutputIterator<char> it_t>
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
        meta::for_each(valid_formats{}, [&] (auto && fmt)
        {
            using fmt_type = remove_cvref_t<decltype(fmt)>;
            using fmt_tag  = typename fmt_type::format_tag;

            for (auto const & ext : fmt_tag::file_extensions)
            {
                if (std::ranges::equal(ext, extension))
                {
                    format = fmt_type{};
                    format_found = true;
                    return;
                }
            }
        });
    }

    if (!format_found)
        throw unhandled_extension_error("No valid format found for this extension.");
}

} // namespace seqan3::detail
