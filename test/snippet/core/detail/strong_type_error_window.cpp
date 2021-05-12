#include <seqan3/core/detail/strong_type.hpp>

struct error : seqan3::detail::strong_type<uint8_t, error>
{
    using seqan3::detail::strong_type<uint8_t, error>::strong_type;
};

struct window_size : seqan3::detail::strong_type<uint8_t, window_size>
{
    using seqan3::detail::strong_type<uint8_t, window_size>::strong_type;
};
