#include <seqan3/alignment/pairwise/alignment_configurator.hpp>

using namespace seqan3;
using namespace seqan3::detail;

int main()
{
    using sequences_t = std::vector<std::pair<std::string, std::string>>;
    using config_t = decltype(align_cfg::edit);

//! [result]
    using first_seq_t = std::tuple_element_t<0, value_type_t<std::remove_reference_t<sequences_t>>>;
    using second_seq_t = std::tuple_element_t<1, value_type_t<std::remove_reference_t<sequences_t>>>;

    // Select the result type based on the sequences and the configuration.
    using result_t = alignment_result<typename align_result_selector<std::remove_reference_t<first_seq_t>,
                                                                 std::remove_reference_t<second_seq_t>,
                                                                 config_t>::type>;
    // Define the function wrapper type.
    using function_wrapper_t = std::function<result_t(first_seq_t &, second_seq_t &)>;
//! [result]

    static_assert(is_type_specialisation_of_v<function_wrapper_t, std::function>);
}
