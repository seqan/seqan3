#include <seqan3/core/detail/strong_type.hpp>

using seqan3::operator|;

struct error : seqan3::detail::strong_type<unsigned, error, seqan3::detail::strong_type_skill::decrement |
                                                            seqan3::detail::strong_type_skill::increment>
{
    using seqan3::detail::strong_type<unsigned, error, seqan3::detail::strong_type_skill::decrement |
                                                       seqan3::detail::strong_type_skill::increment>::strong_type;
};

int main()
{
    error e{4};
    --e;
    ++e;
}
