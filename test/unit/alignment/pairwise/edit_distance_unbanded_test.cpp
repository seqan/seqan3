#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>

#include <seqan3/range/view/to_char.hpp>

#include "fixture/global_edit_distance_unbanded.hpp"
#include "fixture/global_edit_distance_max_errors_unbanded.hpp"
#include "fixture/semi_global_edit_distance_unbanded.hpp"
#include "fixture/semi_global_edit_distance_max_errors_unbanded.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::test::alignment::fixture;

template <typename word_t = uint64_t, typename is_semi_global_t = std::false_type>
struct test_traits_type
{
    using word_type = word_t;
    using is_semi_global_type = is_semi_global_t;
};

template <auto _fixture, typename word_t>
struct param : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }

    using word_type = word_t;
    using is_semi_global_type = std::false_type;
};

template <auto _fixture, typename word_t>
struct param_semi : public param<_fixture, word_t>
{
    using is_semi_global_type = std::true_type;
};

template <typename param_t>
class edit_distance_unbanded : public param_t
{};

TYPED_TEST_CASE_P(edit_distance_unbanded);

using global_edit_distance_unbanded_types
    = ::testing::Types<
        param<&global::edit_distance::unbanded::dna4_01, uint8_t>,
        param<&global::edit_distance::unbanded::dna4_01, uint16_t>,
        param<&global::edit_distance::unbanded::dna4_01, uint32_t>,
        param<&global::edit_distance::unbanded::dna4_01, uint64_t>,

        param<&global::edit_distance::unbanded::dna4_01T, uint8_t>,
        param<&global::edit_distance::unbanded::dna4_01T, uint16_t>,
        param<&global::edit_distance::unbanded::dna4_01T, uint32_t>,
        param<&global::edit_distance::unbanded::dna4_01T, uint64_t>,

        param<&global::edit_distance::unbanded::dna4_02, uint8_t>,
        param<&global::edit_distance::unbanded::dna4_02, uint16_t>,
        param<&global::edit_distance::unbanded::dna4_02, uint32_t>,
        param<&global::edit_distance::unbanded::dna4_02, uint64_t>,

        param<&global::edit_distance::unbanded::aa27_01, uint8_t>,
        param<&global::edit_distance::unbanded::aa27_01, uint16_t>,
        param<&global::edit_distance::unbanded::aa27_01, uint32_t>,
        param<&global::edit_distance::unbanded::aa27_01, uint64_t>,

        param<&global::edit_distance::unbanded::aa27_01T, uint8_t>,
        param<&global::edit_distance::unbanded::aa27_01T, uint16_t>,
        param<&global::edit_distance::unbanded::aa27_01T, uint32_t>,
        param<&global::edit_distance::unbanded::aa27_01T, uint64_t>
    >;

using semi_global_edit_distance_unbanded_types
    = ::testing::Types<
        param_semi<&semi_global::edit_distance::unbanded::dna4_01, uint8_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_01, uint16_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_01, uint32_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_01, uint64_t>,

        param_semi<&semi_global::edit_distance::unbanded::dna4_01T, uint8_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_01T, uint16_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_01T, uint32_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_01T, uint64_t>,

        param_semi<&semi_global::edit_distance::unbanded::dna4_02, uint8_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_02, uint16_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_02, uint32_t>,
        param_semi<&semi_global::edit_distance::unbanded::dna4_02, uint64_t>,

        param_semi<&semi_global::edit_distance::unbanded::aa27_01, uint8_t>,
        param_semi<&semi_global::edit_distance::unbanded::aa27_01, uint16_t>,
        param_semi<&semi_global::edit_distance::unbanded::aa27_01, uint32_t>,
        param_semi<&semi_global::edit_distance::unbanded::aa27_01, uint64_t>,

        param_semi<&semi_global::edit_distance::unbanded::aa27_01T, uint8_t>,
        param_semi<&semi_global::edit_distance::unbanded::aa27_01T, uint16_t>,
        param_semi<&semi_global::edit_distance::unbanded::aa27_01T, uint32_t>,
        param_semi<&semi_global::edit_distance::unbanded::aa27_01T, uint64_t>
    >;

using global_edit_distance_max_errors_unbanded_types
    = ::testing::Types<
        param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint8_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint16_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint32_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint64_t>,

        param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint8_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint16_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint32_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint64_t>,

        param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint8_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint16_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint32_t>,
        param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint64_t>,

        param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint8_t>,
        param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint16_t>,
        param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint32_t>,
        param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint64_t>,

        param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint8_t>,
        param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint16_t>,
        param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint32_t>,
        param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint64_t>
    >;

using semi_global_edit_distance_max_errors_unbanded_types
    = ::testing::Types<
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint8_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint16_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint32_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint64_t>,

        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint8_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint16_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint32_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint64_t>,

        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint8_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint16_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint32_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint64_t>,

        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint8_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint16_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint32_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint64_t>,

        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint8_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint16_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint32_t>,
        param_semi<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint64_t>
    >;

template <typename TypeParam, typename database_t, typename query_t, typename align_cfg_t>
auto edit_distance(database_t && database, query_t && query, align_cfg_t && align_cfg)
{
    using traits_t = test_traits_type<typename TypeParam::word_type, typename TypeParam::is_semi_global_type>;
    using algorithm_t = pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_cfg_t, traits_t>;

    auto result = alignment_result{detail::alignment_result_value_type{}};
    auto alignment = algorithm_t{database, query, align_cfg};

    // compute alignment
    alignment(result);
    return alignment;
}

TYPED_TEST_P(edit_distance_unbanded, score)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam>(database, query, align_cfg);
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded, score_matrix)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam>(database, query, align_cfg);
    auto score_matrix = alignment.score_matrix();

    EXPECT_EQ(score_matrix.cols(), database.size()+1);
    EXPECT_EQ(score_matrix.rows(), query.size()+1);
    EXPECT_EQ(score_matrix, fixture.score_matrix());
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded, trace_matrix)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam>(database, query, align_cfg);
    auto trace_matrix = alignment.trace_matrix();

    EXPECT_EQ(trace_matrix.cols(), database.size()+1);
    EXPECT_EQ(trace_matrix.rows(), query.size()+1);
    EXPECT_EQ(trace_matrix, fixture.trace_matrix());
}

TYPED_TEST_P(edit_distance_unbanded, back_coordinate)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam>(database, query, align_cfg);
    auto back_coordinate = alignment.back_coordinate();

    EXPECT_EQ(back_coordinate, fixture.back_coordinate);
}

TYPED_TEST_P(edit_distance_unbanded, front_coordinate)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam>(database, query, align_cfg);
    auto front_coordinate = alignment.front_coordinate();

    EXPECT_EQ(front_coordinate, fixture.front_coordinate);
}

TYPED_TEST_P(edit_distance_unbanded, alignment)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam>(database, query, align_cfg);

    auto && [gapped_database, gapped_query] = alignment.alignment();
    EXPECT_EQ(std::string{gapped_database | view::to_char}, fixture.aligned_sequence1);
    EXPECT_EQ(std::string{gapped_query | view::to_char}, fixture.aligned_sequence2);
}

REGISTER_TYPED_TEST_CASE_P(edit_distance_unbanded, score, score_matrix, trace_matrix, back_coordinate, front_coordinate, alignment);

// work around a bug that you can't specify more than 50 template arguments to ::testing::types
INSTANTIATE_TYPED_TEST_CASE_P(global, edit_distance_unbanded, global_edit_distance_unbanded_types);
INSTANTIATE_TYPED_TEST_CASE_P(semi_global, edit_distance_unbanded, semi_global_edit_distance_unbanded_types);
INSTANTIATE_TYPED_TEST_CASE_P(global_max_errors, edit_distance_unbanded, global_edit_distance_max_errors_unbanded_types);
INSTANTIATE_TYPED_TEST_CASE_P(semi_global_max_errors, edit_distance_unbanded, semi_global_edit_distance_max_errors_unbanded_types);
