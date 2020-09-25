#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/output.hpp>

int main()
{
    // Only return the reference id where a query matched the reference:
    seqan3::configuration const cfg1 = seqan3::search_cfg::output_reference_id{};

    // Same as the default:
    seqan3::configuration const cfg2 = seqan3::search_cfg::output_query_id{} |
                                       seqan3::search_cfg::output_reference_id{} |
                                       seqan3::search_cfg::output_reference_begin_position{};

    // Only return cursors of the index.
    seqan3::configuration const cfg3 = seqan3::search_cfg::output_index_cursor{};
    return 0;
}
