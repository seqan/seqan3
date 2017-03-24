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

template <typename file_base_traits>
class file_base
{
public:
    /* types */
    using stream_type = typename file_base_traits::stream_type;
    using valid_format_types = typename file_base_traits::valid_format_types;

protected:

    /* constructors */
    // constructor with arg
    file_base(std::experimental::filesystem::path _file_name)
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
    file_base(file_base &&) = default;
    file_base & operator=(file_base &&) = default;

    ~file_base() = default;

    /* member variables */
    stream_type stream;
    valid_format_types format;

    /* member functions */
    std::experimental::filesystem::path select_compression_format(std::experimental::filesystem::path & file_name);
    template <size_t index>
    void select_format(std::experimental::filesystem::path const & ext);
};

// ------------------------------------------------------------------
// protected functions
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

        if (std::find(file_exts.begin(), file_exts.end(), ext) != file_exts.end())
            format = std::variant_alternative_t<index, valid_format_types>{};
        else
            select_format<index+1>(ext);
    }
}

} // namespace detail
