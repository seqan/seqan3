// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \brief Provides various utility functions.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/std/iterator>

namespace seqan3::detail
{

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
                filesystem::path const & file_name)
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

            for (auto const & ext : fmt_type::file_extensions)
            {
                if (ranges::equal(ext, extension))
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
