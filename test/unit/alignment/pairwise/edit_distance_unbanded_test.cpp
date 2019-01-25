
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
using namespace seqan3::fixture;

template <typename word_t = uint64_t>
struct test_traits_type
{
    using word_type = word_t;
};

template <auto _fixture, typename word_t>
struct param : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }

    using word_type = word_t;
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
        param<&semi_global::edit_distance::unbanded::dna4_01, uint8_t>,
        param<&semi_global::edit_distance::unbanded::dna4_01, uint16_t>,
        param<&semi_global::edit_distance::unbanded::dna4_01, uint32_t>,
        param<&semi_global::edit_distance::unbanded::dna4_01, uint64_t>,

        param<&semi_global::edit_distance::unbanded::dna4_01T, uint8_t>,
        param<&semi_global::edit_distance::unbanded::dna4_01T, uint16_t>,
        param<&semi_global::edit_distance::unbanded::dna4_01T, uint32_t>,
        param<&semi_global::edit_distance::unbanded::dna4_01T, uint64_t>,

        param<&semi_global::edit_distance::unbanded::dna4_02, uint8_t>,
        param<&semi_global::edit_distance::unbanded::dna4_02, uint16_t>,
        param<&semi_global::edit_distance::unbanded::dna4_02, uint32_t>,
        param<&semi_global::edit_distance::unbanded::dna4_02, uint64_t>,

        param<&semi_global::edit_distance::unbanded::aa27_01, uint8_t>,
        param<&semi_global::edit_distance::unbanded::aa27_01, uint16_t>,
        param<&semi_global::edit_distance::unbanded::aa27_01, uint32_t>,
        param<&semi_global::edit_distance::unbanded::aa27_01, uint64_t>,

        param<&semi_global::edit_distance::unbanded::aa27_01T, uint8_t>,
        param<&semi_global::edit_distance::unbanded::aa27_01T, uint16_t>,
        param<&semi_global::edit_distance::unbanded::aa27_01T, uint32_t>,
        param<&semi_global::edit_distance::unbanded::aa27_01T, uint64_t>
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
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint8_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint16_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint32_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint64_t>,

        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint8_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint16_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint32_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint64_t>,

        param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint8_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint16_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint32_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint64_t>,

        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint8_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint16_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint32_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint64_t>,

        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint8_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint16_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint32_t>,
        param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint64_t>
    >;

template <typename word_type, typename database_t, typename query_t, typename align_cfg_t>
auto edit_distance(database_t && database, query_t && query, align_cfg_t && align_cfg)
{
    using traits_t = test_traits_type<word_type>;
    using algorithm_t = pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_cfg_t, traits_t>;

    auto result = align_result{detail::align_result_value_type{}};
    auto alignment = algorithm_t{database, query, align_cfg};

    // compute alignment
    alignment(result);
    return alignment;
}

TYPED_TEST_P(edit_distance_unbanded, score)
{
    using word_type = typename TypeParam::word_type;
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<word_type>(database, query, align_cfg);
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded, score_matrix)
{
    using word_type = typename TypeParam::word_type;
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<word_type>(database, query, align_cfg);
    auto score_matrix = alignment.score_matrix();

    EXPECT_EQ(score_matrix.cols(), database.size()+1);
    EXPECT_EQ(score_matrix.rows(), query.size()+1);
    EXPECT_EQ(score_matrix, fixture.score_matrix);
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded, trace_matrix)
{
    using word_type = typename TypeParam::word_type;
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<word_type>(database, query, align_cfg);
    auto trace_matrix = alignment.trace_matrix();
    auto begin_coordinate = alignment.begin_coordinate();
    auto end_coordinate = alignment.end_coordinate();

    EXPECT_EQ(trace_matrix.cols(), database.size()+1);
    EXPECT_EQ(trace_matrix.rows(), query.size()+1);
    EXPECT_EQ(begin_coordinate.first_seq_pos, fixture.begin_coordinate.first_seq_pos);
    EXPECT_EQ(begin_coordinate.second_seq_pos, fixture.begin_coordinate.second_seq_pos);
    EXPECT_EQ(end_coordinate.first_seq_pos, fixture.end_coordinate.first_seq_pos);
    EXPECT_EQ(end_coordinate.second_seq_pos, fixture.end_coordinate.second_seq_pos);
    EXPECT_EQ(trace_matrix, fixture.trace_matrix);
    EXPECT_EQ(alignment.score(), fixture.score);

    auto && [gapped_database, gapped_query] = alignment.alignment();
    EXPECT_EQ(std::string{gapped_database | view::to_char}, fixture.gapped_sequence1);
    EXPECT_EQ(std::string{gapped_query | view::to_char}, fixture.gapped_sequence2);
}

TYPED_TEST_P(edit_distance_unbanded, trace)
{
    using word_type = typename TypeParam::word_type;
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<word_type>(database, query, align_cfg);

    auto && [gapped_database, gapped_query] = alignment.alignment();
    EXPECT_EQ(std::string{gapped_database | view::to_char}, fixture.gapped_sequence1);
    EXPECT_EQ(std::string{gapped_query | view::to_char}, fixture.gapped_sequence2);
}

REGISTER_TYPED_TEST_CASE_P(edit_distance_unbanded, score, score_matrix, trace_matrix, trace);

// work around a bug that you can't specify more than 50 template arguments to ::testing::types
INSTANTIATE_TYPED_TEST_CASE_P(global, edit_distance_unbanded, global_edit_distance_unbanded_types);
INSTANTIATE_TYPED_TEST_CASE_P(semi_global, edit_distance_unbanded, semi_global_edit_distance_unbanded_types);
INSTANTIATE_TYPED_TEST_CASE_P(global_max_errors, edit_distance_unbanded, global_edit_distance_max_errors_unbanded_types);
INSTANTIATE_TYPED_TEST_CASE_P(semi_global_max_errors, edit_distance_unbanded, semi_global_edit_distance_max_errors_unbanded_types);
