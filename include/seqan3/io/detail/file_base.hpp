#include <experimental/filesystem>
#include <tuple>
#include <variant>

using namespace std::experimental::filesystem;

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
    file_base(filesystem::path _file_name)
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
    filesystem::path select_compression_format(filesystem::path & file_name);
    template <size_t index>
    void select_format(filesystem::path const & ext);
};

// ------------------------------------------------------------------
// protected functions
// ------------------------------------------------------------------

template <typename file_base_traits>
filesystem::path file_base<file_base_traits>::select_compression_format(filesystem::path & file_name)
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
inline void file_base<file_base_traits>::select_format(filesystem::path const & ext)
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
