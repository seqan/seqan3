// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================
// Author: Svenja Mehringer <svenja.mehringer@fu-berlin.de>
// ==========================================================================

#include <experimental/filesystem>
#include <tuple>
#include <variant>

namespace seqan3::detail
{

// ==================================================================
// file_base
// ==================================================================

//! The base file format
/*! This class is used to supply all file formats with basic functionality
 * like opening the stream and inferring the compression or file format.
 */
template <typename file_base_traits>
class file_base
{
public:
    /* types */
    using stream_type        = typename file_base_traits::stream_type;        //!< The stream type to write or read from
    using valid_format_types = typename file_base_traits::valid_format_types; //!< The valid format types to choose from

protected:

    //! constructor with file name argument
    /*!
     * Passing a file name (path) as an argument to the constructor will open
     * the stream using this name.
     *
     * Note: The sequence file format will automatically be deduced by
     * the extension file name extension:
     *
     *    -# Check whether a valid compression format was used. Is so, strip the
     *        strip the file name and continue, if not continue with the original
     *        file name
     *
     *    -# Check every valid file format in `valid_format_types` for their
     *        extension identifiers in `file_extensions` and choose the according
     *        format. this function will __throw__ when the format cannot be inferred.
     */
    explicit file_base(std::experimental::filesystem::path _file_name)
    {
        stream.open(_file_name, std::ios::binary); // open stream
        select_format<0>((select_compression_format(_file_name)).extension());
    }
    // TODO add constructor with stream object that needs manual type selection

    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    file_base() = delete;
    file_base(file_base const &) = delete;
    file_base & operator=(file_base const &) = delete;

    // move construction and assignment are defaulted
    file_base(file_base &&) = default; //!< default move constructor
    file_base & operator=(file_base &&) = default; //!< default assignment move constructor
    ~file_base() = default; //!< default deconstructor

    /* member variables */
    stream_type stream; //!< the stream object to read from or write to.
    valid_format_types format; //!< the format object to use for tag dispatching

    /* member functions */
    //! Helper function to select the compression format
    /*!
     * This function goes through the valid_compression_formats of the traits
     * object (see e.g. `sequence_file_in_default_traits`) and compares the
     * extension identifier to infer the compression format.
     * If a compression format is found by a matching extension the file_name
     * is stripped of the extension otherwise it is return directly.
     * \param file_name the file name (path).
     * \return file_name (stripped from valid compression format extension).
     */
    std::experimental::filesystem::path select_compression_format(std::experimental::filesystem::path & file_name);
    //! Helper function to select the file format
    /*!
     * This function iterates over the valid_format_types (std::variant) using the template
     * parameter index.
     * For each type the `file_extensions` list is scanned for the
     * file name extension `ext`. If ext is matched to a format, an instantiation
     * of the type is assigned to the `format` member variable of the file object.
     * If no format can be inferred the function __throws__ an error.
     * \param ext the file name extension (e.g. ".fa").
     */
    template <size_t index>
    void select_format(std::experimental::filesystem::path const & ext);
};

// ==================================================================
// protected functions
// ==================================================================

// ------------------------------------------------------------------
// function select_compression_format
// ------------------------------------------------------------------

template <typename file_base_traits>
std::experimental::filesystem::path file_base<file_base_traits>::select_compression_format(std::experimental::filesystem::path & file_name)
{
    for (auto const & pair : file_base_traits::valid_compression_formats)
    {
        if (file_name.extension() == std::get<0>(pair))
        {
            //TODO:: std::visit([&] (auto & compressor) { stream.push(compressor); }, std::get<1>(pair));
            return file_name.replace_extension(""); // return truncated file name
        }
    }
    return file_name; // return original file name if no compression format was found
}

// ------------------------------------------------------------------
// function select format
// -----------------------------------------------------------------

template <typename file_base_traits>
template <size_t index>
inline void file_base<file_base_traits>::select_format(std::experimental::filesystem::path const & ext)
{
    if constexpr (index == std::variant_size_v<valid_format_types>)
    {
        throw std::runtime_error("No valid format found for this extension");
    }
    else
    {
        auto & file_exts = std::variant_alternative_t<index, valid_format_types>::file_extensions;

        if (std::find(std::begin(file_exts), std::end(file_exts), ext) != std::end(file_exts))
            format = std::variant_alternative_t<index, valid_format_types>{};
        else
            select_format<index+1>(ext);
    }
}

} // namespace detail
