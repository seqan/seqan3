#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/output.hpp>

int main()
{
    // Report hits as positions in the text.
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} |
                                       seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}} |
                                       seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}} |
                                       seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}} |
                                       seqan3::search_cfg::output{seqan3::search_cfg::text_position};

    // Return cursors of the index.
    seqan3::configuration const cfg2 = seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}} |
                                       seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}} |
                                       seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}} |
                                       seqan3::search_cfg::output{seqan3::search_cfg::index_cursor};
    return 0;
}
