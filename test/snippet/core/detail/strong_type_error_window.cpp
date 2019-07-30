#include <seqan3/core/detail/strong_type.hpp>

struct error : seqan3::detail::strong_type<unsigned, error>
{
    using seqan3::detail::strong_type<unsigned, error>::strong_type;
};

struct window_size : seqan3::detail::strong_type<unsigned, window_size>
{
    using seqan3::detail::strong_type<unsigned, window_size>::strong_type;
};
